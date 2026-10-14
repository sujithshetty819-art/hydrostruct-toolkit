"""
Gravity siphon capacity and crest-pressure (cavitation / air-release) check.

For temporary water-supply siphons over a weir or embankment in HDPE pipe.

  * HDPE pipe inner diameter from nominal OD and DR (IPS sizes built in)
  * Darcy-Weisbach friction, Swamee-Jain friction factor
  * Solves the flow for a given static head (upstream WL - outlet level)
  * Pressure at the siphon crest vs vapour pressure, using site atmospheric
    pressure from elevation (ISA). Siphons start losing prime from dissolved
    air release well before vapour pressure, so a margin is applied.

Units: SI. Elevations in m (geodetic), flows m^3/s.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

G = 9.81
RHO = 1000.0
GAMMA = RHO * G

# IPS outside diameters (mm) for common nominal sizes (inches)
IPS_OD_MM = {2: 60.33, 3: 88.90, 4: 114.30, 6: 168.28, 8: 219.08,
             10: 273.05, 12: 323.85, 14: 355.60, 16: 406.40, 18: 457.20,
             20: 508.00, 24: 609.60}


def hdpe_id(nominal_in: int, dr: float) -> float:
    """Average inner diameter (m) of IPS HDPE pipe: ID = OD - 2.12 t (ASTM F714
    average-ID convention), t = OD / DR."""
    od = IPS_OD_MM[nominal_in]
    t = od / dr
    return (od - 2.12 * t) / 1000.0


def water_viscosity(temp_c: float) -> float:
    """Kinematic viscosity of water, m^2/s (Vogel fit, 0-40 C)."""
    mu = 2.414e-5 * 10 ** (247.8 / (temp_c + 273.15 - 140.0))
    return mu / RHO


def vapour_pressure_kpa(temp_c: float) -> float:
    """Saturation vapour pressure of water, kPa (Tetens)."""
    return 0.61078 * math.exp(17.27 * temp_c / (temp_c + 237.3))


def atm_pressure_kpa(elev_m: float) -> float:
    """Standard-atmosphere pressure at elevation, kPa."""
    return 101.325 * (1 - 2.25577e-5 * elev_m) ** 5.25588


def swamee_jain(re: float, rel_rough: float) -> float:
    if re < 2300:
        return 64.0 / re
    return 0.25 / math.log10(rel_rough / 3.7 + 5.74 / re ** 0.9) ** 2


@dataclass
class Siphon:
    D: float                 # inner diameter, m
    L_total: float           # total pipe length inlet -> outlet, m
    L_to_crest: float        # pipe length inlet -> crest, m
    wl_up: float             # upstream water level, m
    z_crest: float           # pipe crown elevation at the crest, m
    z_out: float             # outlet: downstream WL if submerged, else pipe CL
    k_inlet: float = 0.5     # entrance loss
    k_bends_to_crest: float = 0.6
    k_bends_after_crest: float = 0.6
    k_outlet: float = 1.0    # exit loss (submerged outlet = 1.0)
    roughness: float = 0.0015e-3   # HDPE, m
    temp_c: float = 5.0

    @property
    def area(self) -> float:
        return math.pi * self.D ** 2 / 4

    def _f(self, v: float) -> float:
        re = max(v * self.D / water_viscosity(self.temp_c), 1.0)
        return swamee_jain(re, self.roughness / self.D)

    def solve(self) -> dict:
        H = self.wl_up - self.z_out
        if H <= 0:
            raise ValueError("No driving head: outlet is above upstream WL")
        k_minor = (self.k_inlet + self.k_bends_to_crest
                   + self.k_bends_after_crest + self.k_outlet)
        v = 1.0
        for _ in range(100):
            f = self._f(v)
            v_new = math.sqrt(2 * G * H / (f * self.L_total / self.D + k_minor))
            if abs(v_new - v) < 1e-10:
                break
            v = v_new
        f = self._f(v)
        hv = v ** 2 / (2 * G)
        loss_to_crest = (self.k_inlet + self.k_bends_to_crest
                         + f * self.L_to_crest / self.D) * hv

        # Energy: p_crest/gamma = p_atm/gamma + (wl_up - z_crest) - hv - losses
        p_atm = atm_pressure_kpa(self.wl_up) * 1e3
        p_crest_abs = p_atm + GAMMA * (self.wl_up - self.z_crest - hv
                                       - loss_to_crest)
        p_vap = vapour_pressure_kpa(self.temp_c) * 1e3
        head_margin = (p_crest_abs - p_vap) / GAMMA
        max_lift = (p_atm - p_vap) / GAMMA - hv - loss_to_crest

        return {
            "Q_m3s": v * self.area, "V_ms": v, "f": f, "H_static_m": H,
            "p_crest_abs_kPa": p_crest_abs / 1e3,
            "p_crest_gauge_m": (p_crest_abs - p_atm) / GAMMA,
            "margin_above_vapour_m": head_margin,
            "theoretical_max_lift_m": max_lift,
            "actual_lift_m": self.z_crest - self.wl_up,
            # Practical guidance: keep crest absolute head >= ~2.5-3 m above
            # vapour to limit air release / loss of prime.
            "crest_ok": head_margin >= 2.5,
        }

    def flow_vs_upstream_level(self, levels: list[float]) -> list[tuple]:
        rows = []
        for wl in levels:
            s = Siphon(**{**self.__dict__, "wl_up": wl})
            try:
                r = s.solve()
                rows.append((wl, r["Q_m3s"], r["margin_above_vapour_m"],
                             r["crest_ok"]))
            except ValueError:
                rows.append((wl, 0.0, float("nan"), False))
        return rows


if __name__ == "__main__":
    D = hdpe_id(10, 11)
    print(f'10" IPS DR11 HDPE: ID = {D*1000:.1f} mm')
    s = Siphon(D=D, L_total=180.0, L_to_crest=40.0, wl_up=100.0,
               z_crest=104.5, z_out=96.0, temp_c=5.0)
    r = s.solve()
    for k, v in r.items():
        print(f"  {k:26s} {v:.4f}" if isinstance(v, float) else f"  {k:26s} {v}")
    target_q = 0.08
    print(f"\nTarget flow {target_q} m3/s -> "
          f"{'MET' if r['Q_m3s'] >= target_q else 'NOT MET'}")
    print("\nDrawdown sensitivity (upstream WL, Q, crest margin, OK?)")
    for wl, q, m, ok in s.flow_vs_upstream_level([100.0, 99.5, 99.0, 98.5, 98.0]):
        print(f"  {wl:6.1f}  {q:6.3f}  {m:6.2f}  {ok}")
