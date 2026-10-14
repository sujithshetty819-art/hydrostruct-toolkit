"""Hand-check tests. Run with:  python -m unittest discover -s tests"""

import math
import unittest

from hydrostruct.gate_jump_scour import SlideGate
from hydrostruct.rc_section import RCSection
from hydrostruct.rebar_fatigue import (allowable_uniform_range,
                                       cycles_to_failure)
from hydrostruct.siphon import Siphon, atm_pressure_kpa, hdpe_id
from hydrostruct.trashrack import Trashrack


class TestRCSection(unittest.TestCase):
    def setUp(self):
        self.s = RCSection(b=1000, h=450, cover=75, bar_dia=20, spacing=150,
                           fc=35, fy=400)

    def test_flexure_hand_calc(self):
        As = math.pi * 100 * 1000 / 150
        a = 0.85 * As * 400 / (0.7975 * 0.65 * 35 * 1000)
        Mr = 0.85 * As * 400 * (365 - a / 2) / 1e6
        self.assertAlmostEqual(self.s.flexure()["Mr_kNm"], Mr, places=6)

    def test_alpha_beta(self):
        self.assertAlmostEqual(self.s.alpha1, 0.7975)
        self.assertAlmostEqual(self.s.beta1, 0.8825)

    def test_service_stress_equilibrium(self):
        r = self.s.service_stress(65.0)
        C = 0.5 * r["fc_MPa"] * r["k"] * self.s.d * self.s.b
        T = r["fs_MPa"] * self.s.As
        self.assertAlmostEqual(C, T, delta=1e-6 * T)


class TestGate(unittest.TestCase):
    def test_belanger(self):
        r = SlideGate(2.0, 0.5, 3.0, 1.0).analyse()
        y1, fr = r["y1_m"], r["Fr1"]
        self.assertAlmostEqual(
            r["y2_m"], y1 / 2 * (math.sqrt(1 + 8 * fr ** 2) - 1))

    def test_energy_conserved_to_vena_contracta(self):
        g = SlideGate(2.0, 0.5, 3.0, 1.0)
        r = g.analyse()
        E0 = 3.0 + r["q_m2s"] ** 2 / (2 * 9.81 * 9.0)
        E1 = r["y1_m"] + r["v1_ms"] ** 2 / (2 * 9.81)
        self.assertAlmostEqual(E0, E1, places=6)

    def test_high_tailwater_drowns_jump(self):
        self.assertTrue(SlideGate(2.0, 0.3, 3.0, 2.5).analyse()["drowned"])


class TestSiphon(unittest.TestCase):
    def test_hdpe_id(self):
        self.assertAlmostEqual(hdpe_id(12, 17) * 1000, 283.5, delta=0.2)

    def test_sea_level_pressure(self):
        self.assertAlmostEqual(atm_pressure_kpa(0), 101.325)

    def test_energy_balance(self):
        s = Siphon(D=0.2835, L_total=180, L_to_crest=40, wl_up=100,
                   z_crest=104.5, z_out=96)
        r = s.solve()
        hv = r["V_ms"] ** 2 / 19.62
        k = s.k_inlet + s.k_bends_to_crest + s.k_bends_after_crest + s.k_outlet
        loss = (r["f"] * s.L_total / s.D + k) * hv
        self.assertAlmostEqual(loss, 4.0, places=6)


class TestTrashrack(unittest.TestCase):
    def test_max_span_meets_ratio(self):
        r = Trashrack(t=0.01, D=0.075, b=0.04, L=1.0)
        L = r.max_span_for_vibration(1.0, min_ratio=2.0)
        r2 = Trashrack(t=0.01, D=0.075, b=0.04, L=L)
        self.assertAlmostEqual(r2.vibration_check(1.0)["ratio"], 2.0, places=6)


class TestFatigue(unittest.TestCase):
    def test_knee_point(self):
        self.assertAlmostEqual(cycles_to_failure(162.5 / 1.15), 1e6)

    def test_allowable_inverse(self):
        for n in (1e5, 1e8):
            ds = allowable_uniform_range(n)
            self.assertAlmostEqual(cycles_to_failure(ds) / n, 1.0, places=6)


if __name__ == "__main__":
    unittest.main()
