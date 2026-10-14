# hydrostruct-toolkit

These are small, readable Python checks for **hydraulic structures**: walls, slabs, gates, intakes, siphons and pumped-storage works. They cover the questions that come up on real projects and that general-purpose structural packages don't answer.

There are no dependencies. The code uses the Python 3.9+ standard library only. Each module runs on its own with a worked example.

| Module | What it checks | Basis |
|---|---|---|
| `rc_section.py` | RC slab/wall strip: Mr, Vc, ductility, As,min, service steel stress, z-factor, crack width | CSA A23.3-19, Gergely–Lutz |
| `gate_jump_scour.py` | Slide-gate discharge, hydraulic jump free vs drowned, jump length, riprap D50 vs gate opening | Belanger, USBR, Isbash |
| `siphon.py` | HDPE siphon flow, crest pressure vs vapour pressure at site elevation, drawdown sensitivity | Darcy–Weisbach, Swamee–Jain, ISA |
| `trashrack.py` | Rack head loss, bar stress under clogging, vortex-shedding vibration, max spacer spacing | Kirschmer, ASCE intake guidelines |
| `rebar_fatigue.py` | Rebar fatigue under daily reservoir cycling, Miner damage, bent-bar reduction | EN 1992-1-1 §6.8, ACI 215R |

## Quick start

```bash
git clone https://github.com/<you>/hydrostruct-toolkit.git
cd hydrostruct-toolkit
python -m hydrostruct.gate_jump_scour
python -m unittest discover -s tests
```

## Examples

**Will riprap survive downstream of a throttled slide gate?**

```python
from hydrostruct.gate_jump_scour import SlideGate, opening_sweep
gate = SlideGate(width=2.5, opening=0.6, upstream_depth=3.5, tailwater=1.4)
for r in opening_sweep(gate, [0.2, 0.4, 0.8]):
    print(r["opening_m"], r["regime"], round(r["riprap_D50_m"], 2))
```

```
 a (m)      Q   Fr1    y2  TW/y2  V(m/s)  D50(m)  regime
  0.20   2.48  7.45  1.23   1.14    0.71    0.02  drowned jump
  0.30   3.70  6.03  1.47   0.95    8.08    2.73  SWEPT-OUT jump (jet on bed)
```

Once tailwater drops below the sequent depth, the jet sweeps the jump off the apron. Riprap sizing then jumps by two orders of magnitude.

**RC wall base, 1 m strip**

```python
from hydrostruct.rc_section import RCSection, format_report
wall = RCSection(b=1000, h=400, cover=75, bar_dia=20, spacing=150, fc=35, fy=400)
print(format_report("Wall base", wall.check(Mf_kNm=95, Vf_kN=110, Ms_kNm=65)))
```

**Trashrack vibration: how far apart can the spacers be?**

```python
from hydrostruct.trashrack import Trashrack
rack = Trashrack(t=0.010, D=0.075, b=0.040, L=0.45, L_support=1.2)
print(rack.vibration_check(v_approach=0.9, blockage=0.5, min_ratio=2.5))
```

**Pumped-storage rebar fatigue (100-year life)**

```python
from hydrostruct.rebar_fatigue import LoadBlock, miner_damage, bend_reduction
blocks = [LoadBlock("Daily cycle", 75, 2), LoadBlock("Load-follow", 40, 6)]
print(miner_damage(blocks, life_years=100, xi=bend_reduction(4, 1))["D_total"])
```

## Conventions

- SI units throughout. The RC module uses N, mm and MPa, and reports kN and kN·m.
- Every function returns a plain `dict`, so results can go straight into pandas, Excel or a calc report.
- Each check reports a utilisation or an OK/NG flag along with the intermediate values, so a reviewer can follow the calculation.

## Disclaimer

These are screening and verification tools for qualified engineers. They are not a substitute for engineering judgement or project-specific code checks. Verify the clauses and coefficients against the edition of the standard your project uses.

## License

MIT
