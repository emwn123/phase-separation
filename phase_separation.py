"""Periodic 2D Cahn-Hilliard dynamics in dimensionless units."""

import argparse
import csv
import json
import math
from itertools import pairwise
from pathlib import Path

import numpy as np


class CahnHilliard:
    """Stabilized semi-implicit Fourier discretization on a square grid.

    dphi/dt = M laplacian(phi**3 - phi - kappa*laplacian(phi)).
    A linear stabilization S*(phi_new-phi_old) is added to the chemical
    potential. Nonlinear products are evaluated on the grid without dealiasing.
    """

    def __init__(self, n=64, length=64.0, dt=0.5, mobility=1.0,
                 kappa=1.0, stabilization=2.0):
        if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 8 or n % 2:
            raise ValueError("n must be an even integer >= 8")
        for name, value in (("length", length), ("dt", dt), ("mobility", mobility),
                            ("kappa", kappa), ("stabilization", stabilization)):
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        self.n, self.length, self.dt = n, length, dt
        self.mobility, self.kappa, self.stabilization = mobility, kappa, stabilization
        self.dx = length / n
        k = 2 * np.pi * np.fft.fftfreq(n, d=self.dx)
        self.k2 = k[:, None]**2 + k[None, :]**2
        self.denominator = 1 + dt * mobility * self.k2 * (stabilization + kappa * self.k2)

    def _field(self, phi):
        phi = np.asarray(phi, dtype=float)
        if phi.shape != (self.n, self.n) or not np.isfinite(phi).all():
            raise ValueError("phi must be a finite n-by-n array")
        return phi

    def initialize(self, mean=0.0, noise=0.05, seed=7):
        if not math.isfinite(mean) or not math.isfinite(noise) or noise < 0:
            raise ValueError("mean must be finite and noise finite and nonnegative")
        perturbation = np.random.default_rng(seed).uniform(-noise, noise, (self.n, self.n))
        return mean + perturbation - perturbation.mean()

    def step(self, phi):
        phi = self._field(phi)
        with np.errstate(over="raise", invalid="raise"):
            spectrum = np.fft.fft2(phi)
            nonlinear = np.fft.fft2(phi**3 - phi)
            updated = (spectrum + self.dt * self.mobility * self.k2 *
                       (self.stabilization * spectrum - nonlinear)) / self.denominator
            result = np.fft.ifft2(updated).real
        if not np.isfinite(result).all():
            raise FloatingPointError("non-finite field; reduce dt and inspect parameters")
        return result

    def diagnostics(self, phi):
        phi = self._field(phi)
        spectrum = np.fft.fft2(phi)
        # Parseval normalization: mean(|grad phi|^2)=sum(k^2 |phi_hat|^2)/n^4.
        gradient_mean = np.sum(self.k2 * np.abs(spectrum)**2) / self.n**4
        energy = self.length**2 * (np.mean((phi**2 - 1)**2 / 4) +
                                  0.5 * self.kappa * gradient_mean)
        return {"mean": float(phi.mean()), "mass": float(phi.mean() * self.length**2),
                "energy": float(energy), "minimum": float(phi.min()),
                "maximum": float(phi.max()), "variance": float(phi.var())}


def simulate(output, *, n=64, steps=2000, dt=0.5, length=64.0, mobility=1.0,
             kappa=1.0, stabilization=2.0, mean=0.0, noise=0.05, seed=7, every=200):
    """Run a reproducible case and save snapshots, diagnostics and a figure."""
    if not isinstance(steps, int) or steps < 1 or not isinstance(every, int) or every < 1:
        raise ValueError("steps and every must be positive integers")
    output = Path(output)
    # Keep prior runs intact instead of silently overwriting their provenance.
    if output.exists() and any(output.iterdir()):
        raise ValueError("output directory is not empty; choose a new output path")
    output.mkdir(parents=True, exist_ok=True)
    model = CahnHilliard(n, length, dt, mobility, kappa, stabilization)
    phi = model.initialize(mean, noise, seed)
    rows, fields, times = [], [], []
    initial_mean = float(phi.mean())
    for iteration in range(steps + 1):
        if iteration % every == 0 or iteration == steps:
            row = {"step": iteration, "time": iteration * dt, **model.diagnostics(phi)}
            rows.append(row)
            fields.append(phi.copy())
            times.append(iteration * dt)
        if iteration < steps:
            phi = model.step(phi)
    np.savez_compressed(output / "fields.npz", phi=np.array(fields), time=np.array(times),
                        length=length, dx=model.dx)
    with (output / "diagnostics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    configuration = {"n": n, "steps": steps, "dt": dt, "length": length,
                     "mobility": mobility, "kappa": kappa, "stabilization": stabilization,
                     "mean": mean, "noise": noise, "seed": seed, "every": every,
                     "numpy_version": np.__version__}
    (output / "config.json").write_text(json.dumps(configuration, indent=2) + "\n",
                                       encoding="utf-8")
    report = {"max_sampled_mean_drift": max(abs(r["mean"] - initial_mean) for r in rows),
              "initial_energy": rows[0]["energy"], "final_energy": rows[-1]["energy"],
              "sampled_energy_increases": sum(b["energy"] > a["energy"] + 1e-10
                                               for a, b in pairwise(rows)),
              "initial_variance": rows[0]["variance"], "final_variance": rows[-1]["variance"]}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    plot(output, fields, times, rows, length)
    return report


def plot(output, fields, times, rows, length):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selected = np.unique(np.linspace(0, len(fields) - 1, min(4, len(fields)), dtype=int))
    fig, axes = plt.subplots(1, len(selected), figsize=(3.4 * len(selected), 3.2),
                             squeeze=False, constrained_layout=True)
    for ax, index in zip(axes[0], selected):
        im = ax.imshow(fields[index], origin="lower", extent=(0, length, 0, length),
                       cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set(title=f"t = {times[index]:g}", xlabel="x", ylabel="y")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="Order parameter, phi", shrink=0.85)
    fig.savefig(output / "phase-separation.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(8, 3), constrained_layout=True)
    t = [r["time"] for r in rows]
    axes[0].plot(t, [r["energy"] for r in rows])
    axes[0].set(xlabel="Time", ylabel="Free energy")
    axes[1].plot(t, [r["mean"] - rows[0]["mean"] for r in rows])
    axes[1].set(xlabel="Time", ylabel="Mean composition drift")
    fig.savefig(output / "diagnostics.png", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("run"))
    parser.add_argument("--n", type=int, default=64)
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--every", type=int, default=200)
    parser.add_argument("--seed", type=int, default=7)
    for name, default in (("dt", 0.5), ("length", 64.0), ("mobility", 1.0),
                          ("kappa", 1.0), ("stabilization", 2.0), ("mean", 0.0),
                          ("noise", 0.05)):
        parser.add_argument(f"--{name}", type=float, default=default)
    args = vars(parser.parse_args())
    try:
        print(json.dumps(simulate(**args), indent=2))
    except (ValueError, OSError, FloatingPointError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
