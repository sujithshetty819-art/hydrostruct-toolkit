"""
Reinforced concrete rectangular section checks to CSA A23.3-19.

Written for the one-way slabs and cantilever walls typical of hydraulic
structures (fishway walls, flume walls, intake slabs), where serviceability
(steel stress and crack width) usually governs over ultimate strength.

Covers:
  * Factored flexural resistance Mr (equivalent rectangular stress block)
  * Ductility check (c/d limit, Cl. 10.5.2)
  * Minimum flexural reinforcement (Cl. 10.5.1.2)
  * Concrete shear resistance Vc, simplified method (Cl. 11.3.4 / 11.3.6.2)
  * Cracked-elastic service stress in steel
  * Crack-control parameter z (Cl. 10.6.1) and Gergely-Lutz crack width

Units: N, mm, MPa (N/mm^2). Moments in kN*m, shears in kN at the interface.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

PHI_C = 0.65   # concrete resistance factor
PHI_S = 0.85   # reinforcing steel resistance factor
ES = 200_000.0  # MPa


def bar_area(diameter_mm: float) -> float:
    """Nominal area of a single bar, mm^2."""
    return math.pi * diameter_mm ** 2 / 4.0


def as_per_metre(bar_dia_mm: float, spacing_mm: float) -> float:
    """Steel area per metre width for bars at a given spacing, mm^2/m."""
    return bar_area(bar_dia_mm) * 1000.0 / spacing_mm


@dataclass
class RCSection:
    """Singly reinforced rectangular section (per-metre strip for slabs/walls).

    b       : width, mm (use 1000 for a 1 m strip)
    h       : overall thickness, mm
    cover   : clear cover to the tension bars, mm
    bar_dia : tension bar diameter, mm
    spacing : tension bar spacing, mm
    fc      : specified concrete strength f'c, MPa
    fy      : specified yield strength, MPa
    lam     : lambda, density factor (1.0 normal density)
    exposure: 'interior' or 'exterior' (crack-control z limit)
    """

    b: float
    h: float
    cover: float
    bar_dia: float
    spacing: float
    fc: float = 35.0
    fy: float = 400.0
    lam: float = 1.0
    exposure: str = "exterior"
    stirrups: bool = False
    max_agg_mm: float = 20.0

    # ---------- geometry ----------
    @property
    def d(self) -> float:
        return self.h - self.cover - self.bar_dia / 2.0

    @property
    def As(self) -> float:
        return bar_area(self.bar_dia) * self.b / self.spacing

    @property
    def dc(self) -> float:
        """Concrete cover to centroid of the bar closest to the tension face."""
        return self.cover + self.bar_dia / 2.0

    # ---------- material ----------
    @property
    def alpha1(self) -> float:
        return max(0.85 - 0.0015 * self.fc, 0.67)

    @property
    def beta1(self) -> float:
        return max(0.97 - 0.0025 * self.fc, 0.67)

    @property
    def Ec(self) -> float:
        """Modulus of elasticity, Cl. 8.6.2.3 simplified (normal density)."""
        return 4500.0 * math.sqrt(self.fc)

    @property
    def n(self) -> float:
        return ES / self.Ec

    # ---------- ULS flexure ----------
    def flexure(self) -> dict:
        a = PHI_S * self.As * self.fy / (self.alpha1 * PHI_C * self.fc * self.b)
        c = a / self.beta1
        Mr = PHI_S * self.As * self.fy * (self.d - a / 2.0) / 1e6  # kN*m
        c_d_limit = 700.0 / (700.0 + self.fy)
        fcr = 0.6 * self.lam * math.sqrt(self.fc)
        As_min = 0.2 * math.sqrt(self.fc) / self.fy * self.b * self.h
        return {
            "d_mm": self.d,
            "As_mm2": self.As,
            "a_mm": a,
            "c_over_d": c / self.d,
            "c_over_d_limit": c_d_limit,
            "ductile": c / self.d <= c_d_limit,
            "Mr_kNm": Mr,
            "Mcr_kNm": fcr * self.b * self.h ** 2 / 6.0 / 1e6,
            "As_min_mm2": As_min,
            "As_min_ok": self.As >= As_min,
        }

    # ---------- ULS shear (simplified method) ----------
    def shear(self) -> dict:
        dv = max(0.9 * self.d, 0.72 * self.h)
        if self.stirrups:
            beta = 0.18
        elif self.h <= 350.0:
            beta = 0.21
        else:
            sze = 35.0 * dv / (15.0 + self.max_agg_mm) if self.max_agg_mm < 20 else dv
            beta = 230.0 / (1000.0 + sze)
        root_fc = min(math.sqrt(self.fc), 8.0)
        Vc = PHI_C * self.lam * beta * root_fc * self.b * dv / 1e3  # kN
        return {"dv_mm": dv, "beta": beta, "Vc_kN": Vc}

    # ---------- SLS ----------
    def service_stress(self, Ms_kNm: float) -> dict:
        """Cracked transformed section stresses under service moment."""
        rho = self.As / (self.b * self.d)
        rn = rho * self.n
        k = math.sqrt(2.0 * rn + rn ** 2) - rn
        j = 1.0 - k / 3.0
        M = Ms_kNm * 1e6
        fs = M / (self.As * j * self.d)
        fc = 2.0 * M / (k * j * self.b * self.d ** 2)
        return {"k": k, "j": j, "fs_MPa": fs, "fc_MPa": fc,
                "fs_over_fy": fs / self.fy}

    def crack_control(self, Ms_kNm: float) -> dict:
        """CSA z-factor and Gergely-Lutz crack width estimate."""
        fs = self.service_stress(Ms_kNm)["fs_MPa"]
        A = 2.0 * self.dc * self.spacing  # effective tension area per bar
        z = fs * (self.dc * A) ** (1.0 / 3.0)
        z_lim = 30_000.0 if self.exposure == "interior" else 25_000.0
        # Gergely-Lutz: w = 11e-6 * beta_h * fs * (dc*A)^(1/3), fs in MPa -> mm
        # (constant converted from the imperial form: 0.076e-3 in, ksi)
        k = self.service_stress(Ms_kNm)["k"]
        c = k * self.d
        beta_h = (self.h - c) / (self.d - c)
        w = 11.0e-6 * beta_h * fs * (self.dc * A) ** (1.0 / 3.0)
        return {"fs_MPa": fs, "z_N_per_mm": z, "z_limit": z_lim,
                "z_ok": z <= z_lim, "w_mm": w, "beta_h": beta_h}

    # ---------- one-shot report ----------
    def check(self, Mf_kNm: float, Vf_kN: float, Ms_kNm: float,
              w_limit_mm: float = 0.25) -> dict:
        fl, sh, cc = self.flexure(), self.shear(), self.crack_control(Ms_kNm)
        return {
            "flexure_DCR": Mf_kNm / fl["Mr_kNm"],
            "shear_DCR": Vf_kN / sh["Vc_kN"],
            "crack_width_mm": cc["w_mm"],
            "crack_width_ok": cc["w_mm"] <= w_limit_mm,
            "z_ok": cc["z_ok"],
            "ductile": fl["ductile"],
            "As_min_ok": fl["As_min_ok"],
            "detail": {"flexure": fl, "shear": sh, "crack": cc},
        }


def format_report(title: str, res: dict) -> str:
    fl, sh, cc = (res["detail"][k] for k in ("flexure", "shear", "crack"))
    ok = lambda b: "OK" if b else "NG"  # noqa: E731
    lines = [
        f"== {title} ==",
        f"  d = {fl['d_mm']:.0f} mm, As = {fl['As_mm2']:.0f} mm2",
        f"  Mr = {fl['Mr_kNm']:.1f} kN*m   DCR = {res['flexure_DCR']:.2f}  "
        f"{ok(res['flexure_DCR'] <= 1)}",
        f"  Vc = {sh['Vc_kN']:.1f} kN (beta={sh['beta']:.3f})  "
        f"DCR = {res['shear_DCR']:.2f}  {ok(res['shear_DCR'] <= 1)}",
        f"  c/d = {fl['c_over_d']:.3f} <= {fl['c_over_d_limit']:.3f}  "
        f"{ok(res['ductile'])}",
        f"  As,min = {fl['As_min_mm2']:.0f} mm2  {ok(res['As_min_ok'])}",
        f"  fs(service) = {cc['fs_MPa']:.0f} MPa, z = {cc['z_N_per_mm']:.0f} N/mm "
        f"<= {cc['z_limit']:.0f}  {ok(res['z_ok'])}",
        f"  crack width = {cc['w_mm']:.3f} mm  {ok(res['crack_width_ok'])}",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    # Example: 400 mm cantilever channel wall, 20 mm bars @ 150, 75 mm cover
    wall = RCSection(b=1000, h=400, cover=75, bar_dia=20, spacing=150,
                     fc=35, fy=400, exposure="exterior")
    print(format_report("Cantilever wall base (1 m strip)",
                        wall.check(Mf_kNm=95.0, Vf_kN=110.0, Ms_kNm=65.0,
                                   w_limit_mm=0.25)))
