import unittest

import numpy as np

from phase_separation import CahnHilliard


class ModelTests(unittest.TestCase):
    def test_uniform_state(self):
        model = CahnHilliard(n=32)
        phi = np.full((32, 32), 0.3)
        np.testing.assert_allclose(model.step(phi), phi, atol=1e-14)

    def test_mass_and_energy(self):
        model = CahnHilliard(n=32, length=32)
        phi = model.initialize(mean=0.1)
        mean = phi.mean()
        previous = model.diagnostics(phi)["energy"]
        for _ in range(500):
            phi = model.step(phi)
            current = model.diagnostics(phi)["energy"]
            self.assertLessEqual(current, previous + 1e-10)
            self.assertAlmostEqual(phi.mean(), mean, places=12)
            previous = current

    def test_linear_mode_amplification(self):
        # For tiny amplitude about phi=0, growth rate is M*k^2*(1-kappa*k^2).
        model = CahnHilliard(n=32, length=32, dt=0.01)
        k = 2 * np.pi / model.length
        mode = np.cos(k * np.arange(32))
        phi = np.broadcast_to(1e-7 * mode, (32, 32)).copy()
        factor = (1 + model.dt * model.mobility * k**2 * (model.stabilization + 1)) / (
            1 + model.dt * model.mobility * k**2 * (model.stabilization + model.kappa*k**2))
        np.testing.assert_allclose(model.step(phi), phi * factor, atol=1e-19, rtol=1e-11)

    def test_uniform_energy(self):
        model = CahnHilliard(n=16, length=10)
        self.assertAlmostEqual(model.diagnostics(np.zeros((16, 16)))["energy"], 25)
        self.assertAlmostEqual(model.diagnostics(np.ones((16, 16)))["energy"], 0)

    def test_seed_and_mean(self):
        model = CahnHilliard(n=16)
        np.testing.assert_array_equal(model.initialize(seed=1), model.initialize(seed=1))
        self.assertAlmostEqual(model.initialize(mean=0.2).mean(), 0.2, places=14)

    def test_invalid_inputs(self):
        for kwargs in ({"n": 9}, {"dt": 0}, {"kappa": -1}, {"length": np.nan}):
            with self.assertRaises(ValueError):
                CahnHilliard(**kwargs)
        with self.assertRaises(ValueError):
            CahnHilliard().step(np.zeros((8, 8)))


if __name__ == "__main__":
    unittest.main()
