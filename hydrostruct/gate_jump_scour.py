"""
Downstream hydraulics of an underflow (slide/sluice) gate: free vs drowned
hydraulic jump, and riprap stability against the jet.

Typical question this answers: "Will riprap downstream of a partially-open
slide gate survive, or will the jet sweep out the jump and scour the bed?"

Method
  1. Free gate discharge: q = Cc*a*sqrt(2g(H - Cc*a)), solved with the approach
     velocity head included (fixed-point iteration).
  2. Vena contracta depth y1 = Cc*a, Froude number Fr1.
  3. Sequent depth y2 (Belanger).  Compare with tailwater:
       TW >= y2 -> jump is held / drowned at the gate (good for riprap)
       TW <  y2 -> jump is swept downstream; supercritical jet on the bed
  4. Jump length L ~ 6.1*y2 (USBR, valid ~4.5 < Fr1 < 13).
  5. Riprap stone size by Isbash: D = V^2 / (2 g C^2 (Gs - 1)),
     C = 0.86 (high turbulence / jet) or 1.20 (low turbulence).

Units: SI (m, s, m^3/s).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

G = 9.81


@dataclass
class SlideGate:
    width: float            # clear gate width, m
    opening: float          # gate opening a, m
    upstream_depth: float   # H, depth above gate sill, m
    tailwater: float        # downstream depth above sill/apron, m
    cc: float = 0.61        # contraction coefficient (sharp-edged gate)

    def discharge(self) -> float:
        """Free-flow discharge Q, m^3/s (includes approach velocity head)."""
        y1 = self.cc * self.opening
        H = self.upstream_depth
        q = y1 * math.sqrt(2 * G * max(H - y1, 1e-9))
        for _ in range(50):
            E0 = H + q ** 2 / (2 * G * H ** 2)
            q_new = y1 * math.sqrt(2 * G * max(E0 - y1, 1e-9))
            if abs(q_new - q) < 1e-9:
                break
            q = q_new
        return q * self.width

    def analyse(self, gs: float = 2.65) -> dict:
        Q = self.discharge()
        q = Q / self.width
        y1 = self.cc * self.opening
        v1 = q / y1
        fr1 = v1 / math.sqrt(G * y1)
        y2 = y1 / 2.0 * (math.sqrt(1 + 8 * fr1 ** 2) - 1)
        v2 = q / y2
        headloss = (y2 - y1) ** 3 / (4 * y1 * y2)
        drowned = self.tailwater >= y2
        submergence = self.tailwater / y2

        # Velocity attacking the bed: the jet if swept out, else the outgoing
        # flow under a turbulent roller (use v at tailwater, C for turbulent).
        if drowned:
            v_design, C, regime = q / self.tailwater, 0.86, "drowned jump"
        else:
            v_design, C, regime = v1, 0.86, "SWEPT-OUT jump (jet on bed)"
        # Isbash nominal stone diameter, taken as D50.
        # Drowned case is a screening value only: free-flow gate equation is
        # not valid once the jump submerges the gate (check submerged Cd).
        d50 = v_design ** 2 / (2 * G * C ** 2 * (gs - 1))

        return {
            "Q_m3s": Q, "q_m2s": q, "y1_m": y1, "v1_ms": v1, "Fr1": fr1,
            "y2_m": y2, "v2_ms": v2, "tailwater_m": self.tailwater,
            "TW_over_y2": submergence, "drowned": drowned, "regime": regime,
            "jump_length_m": 6.1 * y2, "headloss_m": headloss,
            "v_design_ms": v_design, "riprap_D50_m": d50,
            "note": _froude_note(fr1),
        }


def _froude_note(fr: float) -> str:
    if fr < 1.7:
        return "undular jump"
    if fr < 2.5:
        return "weak jump"
    if fr < 4.5:
        return "oscillating jump - wave action; avoid for riprap aprons"
    if fr < 9.0:
        return "steady jump - best basin performance"
    return "strong jump - high bed attack, check concrete apron"


def opening_sweep(gate: SlideGate, openings: list[float]) -> list[dict]:
    """Run the analysis over a range of gate openings (throttling study)."""
    out = []
    for a in openings:
        g = SlideGate(gate.width, a, gate.upstream_depth, gate.tailwater, gate.cc)
        r = g.analyse()
        r["opening_m"] = a
        out.append(r)
    return out


if __name__ == "__main__":
    # 2.5 m wide slide gate, 3.5 m head, 1.4 m tailwater
    gate = SlideGate(width=2.5, opening=0.6, upstream_depth=3.5, tailwater=1.4)
    print(f"{'a (m)':>6} {'Q':>6} {'Fr1':>5} {'y2':>5} {'TW/y2':>6} "
          f"{'V(m/s)':>7} {'D50(m)':>7}  regime")
    for r in opening_sweep(gate, [0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.4]):
        print(f"{r['opening_m']:6.2f} {r['Q_m3s']:6.2f} {r['Fr1']:5.2f} "
              f"{r['y2_m']:5.2f} {r['TW_over_y2']:6.2f} {r['v_design_ms']:7.2f} "
              f"{r['riprap_D50_m']:7.2f}  {r['regime']}")
