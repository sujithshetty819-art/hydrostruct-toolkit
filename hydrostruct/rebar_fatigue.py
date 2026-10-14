"""
Reinforcing-steel fatigue under cyclic hydraulic loading (e.g. pumped-hydro
storage reservoirs and intakes, where water level cycles 1-4 times per day
for a 100-year design life).

Method: Eurocode 2 (EN 1992-1-1, Cl. 6.8.4) S-N curves with Palmgren-Miner
cumulative damage, plus the ACI 215R stress-range screening limit.

  N(dsigma) = N* (dsigma_Rsk / (gamma_s,fat * gamma_F,fat * dsigma))^k
  k = k1 above the knee point, k2 below.

Detail categories (Table 6.3N, recommended values):
  straight_bent : N*=1e6, k1=5, k2=9, dsigma_Rsk=162.5 MPa
  welded        : N*=1e7, k1=3, k2=5, dsigma_Rsk=58.5 MPa
  coupler       : N*=1e7, k1=3, k2=5, dsigma_Rsk=35.0 MPa
Bent bars: dsigma_Rsk reduced by xi = 0.35 + 0.026 * D_mandrel / phi.
"""

from __future__ import annotations

from dataclasses import dataclass

SN_CURVES = {
    "straight_bent": (1e6, 5, 9, 162.5),
    "welded": (1e7, 3, 5, 58.5),
    "coupler": (1e7, 3, 5, 35.0),
}
GAMMA_S_FAT = 1.15
GAMMA_F_FAT = 1.0


def bend_reduction(mandrel_dia: float, bar_dia: float) -> float:
    return min(1.0, 0.35 + 0.026 * mandrel_dia / bar_dia)


def cycles_to_failure(dsigma: float, detail: str = "straight_bent",
                      xi: float = 1.0) -> float:
    n_star, k1, k2, d_rsk = SN_CURVES[detail]
    resist = xi * d_rsk / GAMMA_S_FAT
    act = GAMMA_F_FAT * dsigma
    if act <= 0:
        return float("inf")
    k = k1 if act >= resist else k2
    return n_star * (resist / act) ** k


def aci215_limit(fmin_mpa: float) -> float:
    """ACI 215R straight-bar stress range limit, MPa (fmin + tension)."""
    return 161.0 - 0.33 * fmin_mpa


@dataclass
class LoadBlock:
    name: str
    dsigma_mpa: float        # steel stress range per cycle (service, cracked)
    cycles_per_day: float
    fmin_mpa: float = 0.0


def miner_damage(blocks: list[LoadBlock], life_years: float,
                 detail: str = "straight_bent", xi: float = 1.0) -> dict:
    rows, total = [], 0.0
    for blk in blocks:
        n = blk.cycles_per_day * 365.25 * life_years
        N = cycles_to_failure(blk.dsigma_mpa, detail, xi)
        d = n / N
        total += d
        rows.append({"block": blk.name, "dsigma": blk.dsigma_mpa, "n": n,
                     "N": N, "damage": d,
                     "aci215_ok": blk.dsigma_mpa <= aci215_limit(blk.fmin_mpa)})
    return {"rows": rows, "D_total": total, "ok": total <= 1.0,
            "life_reserve": (1 / total) if total > 0 else float("inf")}


def allowable_uniform_range(n_cycles: float, detail: str = "straight_bent",
                            xi: float = 1.0) -> float:
    """Stress range (MPa) that gives D = 1.0 for n_cycles of constant amplitude."""
    n_star, k1, k2, d_rsk = SN_CURVES[detail]
    resist = xi * d_rsk / GAMMA_S_FAT
    k = k1 if n_cycles <= n_star else k2
    return resist * (n_star / n_cycles) ** (1 / k) / GAMMA_F_FAT


if __name__ == "__main__":
    spectrum = [
        LoadBlock("Full generate/pump cycle", 75.0, 2.0, fmin_mpa=40),
        LoadBlock("Partial load-follow cycle", 40.0, 6.0, fmin_mpa=60),
        LoadBlock("Wave/surge in headpond", 12.0, 2000.0, fmin_mpa=80),
    ]
    for label, xi in [("straight bars", 1.0),
                      ("bent bars, 4phi mandrel", bend_reduction(4, 1))]:
        res = miner_damage(spectrum, life_years=100, xi=xi)
        print(f"\n--- {label} (xi = {xi:.3f}) ---")
        print(f"{'block':28s} {'dsig':>6} {'n':>10} {'N':>10} {'D':>8}  ACI215")
        for r in res["rows"]:
            print(f"{r['block']:28s} {r['dsigma']:6.1f} {r['n']:10.3g} "
                  f"{r['N']:10.3g} {r['damage']:8.4f}  "
                  f"{'OK' if r['aci215_ok'] else 'NG'}")
        print(f"Total Miner damage D = {res['D_total']:.3f} -> "
              f"{'OK' if res['ok'] else 'NG'}")
    n100 = 2 * 365.25 * 100
    print(f"\nAllowable constant-amplitude range for {n100:,.0f} cycles: "
          f"{allowable_uniform_range(n100):.0f} MPa (straight bars)")
