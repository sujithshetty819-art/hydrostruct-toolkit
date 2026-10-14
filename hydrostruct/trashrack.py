"""
Trashrack design: head loss, bar strength under clogging, and flow-induced
vibration (vortex-shedding lock-in) screening.

Bar geometry convention (looking down the flow):
    t = bar thickness (the face the flow hits), m
    D = bar depth in the flow direction, m
    b = clear spacing between bars, m
    L = span of bar between lateral supports, m

References
  * Kirschmer head loss: h = kf (t/b)^(4/3) sin(theta) V^2/2g
  * ASCE (1995) Guidelines for Design of Intakes for Hydroelectric Plants:
    keep bar natural frequency in water >= 1.5x the vortex-shedding frequency
    (a larger margin is common for reversing-flow pumped-storage intakes).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

G = 9.81
RHO_W = 1000.0
RHO_STEEL = 7850.0
E_STEEL = 200e9

KIRSCHMER_KF = {
    "rectangular": 2.42,
    "rounded_upstream": 1.83,
    "circular": 1.79,
    "rounded_both": 1.67,
    "airfoil": 0.76,
}

# Mode-1 frequency coefficients  f = (C / 2pi) sqrt(EI / (m L^4))
END_FIXITY_C = {"pinned": math.pi ** 2, "fixed": 22.37, "fixed_pinned": 15.42}


@dataclass
class Trashrack:
    t: float
    D: float
    b: float
    L: float                        # span between lateral spacers (vibration), m
    L_support: float = 0.0          # span between support beams (bending); 0 -> L
    theta_deg: float = 90.0         # rack inclination from horizontal
    shape: str = "rectangular"
    fixity: str = "pinned"
    fy: float = 300e6               # Pa
    strouhal: float = 0.22          # based on thickness t

    @property
    def open_ratio(self) -> float:
        return self.b / (self.b + self.t)

    def head_loss(self, v_approach: float) -> float:
        kf = KIRSCHMER_KF[self.shape]
        return (kf * (self.t / self.b) ** (4 / 3)
                * math.sin(math.radians(self.theta_deg))
                * v_approach ** 2 / (2 * G))

    def net_velocity(self, v_approach: float, blockage: float = 0.0) -> float:
        return v_approach / (self.open_ratio * (1.0 - blockage))

    # ----- strength under differential head (clogged rack) -----
    def bar_stress(self, dh: float) -> dict:
        w = RHO_W * G * dh * (self.t + self.b)   # N/m on one bar
        coeff = {"pinned": 1 / 8, "fixed": 1 / 12, "fixed_pinned": 1 / 8}
        span = self.L_support or self.L
        M = coeff[self.fixity] * w * span ** 2
        S = self.t * self.D ** 2 / 6
        sigma = M / S
        return {"w_N_per_m": w, "M_Nm": M, "sigma_MPa": sigma / 1e6,
                "utilisation": sigma / (0.9 * self.fy)}

    # ----- vibration -----
    def natural_frequency(self) -> dict:
        """Weak-axis (across t) mode-1 frequency in air and in water."""
        I = self.D * self.t ** 3 / 12
        m_s = RHO_STEEL * self.D * self.t
        # added mass: plate of width D moving normal to its face
        m_a = RHO_W * math.pi * self.D ** 2 / 4
        C = END_FIXITY_C[self.fixity]
        f_air = C / (2 * math.pi) * math.sqrt(E_STEEL * I / (m_s * self.L ** 4))
        f_wat = C / (2 * math.pi) * math.sqrt(E_STEEL * I /
                                               ((m_s + m_a) * self.L ** 4))
        return {"f_air_Hz": f_air, "f_water_Hz": f_wat,
                "added_mass_ratio": m_a / m_s}

    def vibration_check(self, v_approach: float, blockage: float = 0.0,
                        min_ratio: float = 1.5) -> dict:
        vn = self.net_velocity(v_approach, blockage)
        fs = self.strouhal * vn / self.t
        fn = self.natural_frequency()["f_water_Hz"]
        L_max = self.max_span_for_vibration(v_approach, blockage, min_ratio)
        return {"V_net_ms": vn, "f_shed_Hz": fs, "f_n_water_Hz": fn,
                "ratio": fn / fs, "ok": fn / fs >= min_ratio,
                "max_unbraced_span_m": L_max}

    def max_span_for_vibration(self, v_approach, blockage=0.0, min_ratio=1.5):
        """Longest span between lateral spacers that keeps fn/fs >= min_ratio."""
        vn = self.net_velocity(v_approach, blockage)
        fs = self.strouhal * vn / self.t
        f_target = min_ratio * fs
        fn_L = self.natural_frequency()["f_water_Hz"] * self.L ** 2  # fn*L^2 const
        return math.sqrt(fn_L / f_target)


if __name__ == "__main__":
    # 10 x 75 mm flat bars, 40 mm clear, lateral spacers at 450 mm
    rack = Trashrack(t=0.010, D=0.075, b=0.040, L=0.45, L_support=1.2,
                     theta_deg=75,
                     shape="rectangular", fixity="pinned", fy=300e6)
    V = 0.9
    print(f"Open ratio          : {rack.open_ratio:.3f}")
    print(f"Head loss (clean)   : {rack.head_loss(V)*1000:.1f} mm")
    s = rack.bar_stress(dh=2.0)
    print(f"Bar stress @2.0 m dh: {s['sigma_MPa']:.0f} MPa "
          f"(util {s['utilisation']:.2f})")
    for blk in (0.0, 0.25, 0.5):
        v = rack.vibration_check(V, blockage=blk, min_ratio=2.5)
        print(f"Blockage {blk:4.0%}: Vnet={v['V_net_ms']:.2f} m/s  "
              f"fs={v['f_shed_Hz']:.1f} Hz  fn={v['f_n_water_Hz']:.1f} Hz  "
              f"ratio={v['ratio']:.2f} {'OK' if v['ok'] else 'NG'}  "
              f"Lmax={v['max_unbraced_span_m']:.2f} m")
