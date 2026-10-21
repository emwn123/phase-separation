"""Reproduce composition cases and a matched-time time-step comparison."""

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

from phase_separation import simulate

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def run_cases(output):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("choose an empty output directory")
    output.mkdir(parents=True, exist_ok=True)
    cases = [("symmetric", 0.0, 2000), ("off-critical", 0.3, 2000),
             ("stable-mixture", 0.8, 1000)]
    summary = {}
    final_fields = []
    for name, mean, steps in cases:
        summary[name] = simulate(output / name, mean=mean, steps=steps, every=100)
        with np.load(output / name / "fields.npz") as data:
            final_fields.append(data["phi"][-1].copy())

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.5), constrained_layout=True)
    for ax, (name, mean, steps), phi in zip(axes, cases, final_fields):
        im = ax.imshow(phi, origin="lower", extent=(0, 64, 0, 64), cmap="RdBu_r",
                       vmin=-1, vmax=1)
        ax.set(title=f"{name}\nmean = {mean:g}, t = {steps * 0.5:g}", xlabel="x", ylabel="y")
    fig.colorbar(im, ax=axes.tolist(), label="Order parameter, phi", shrink=0.8)
    fig.savefig(output / "composition-cases.png", dpi=180)
    plt.close(fig)

    # Same initial field, domain, grid, stabilization and final time for each dt.
    # The smallest dt is a numerical reference, not an exact solution.
    dt_fields, dt_energy = {}, {}
    for dt in (0.5, 0.25, 0.125, 0.0625):
        name = f"dt-{dt:g}"
        summary[name] = simulate(output / name, n=32, length=32, dt=dt,
                                 steps=round(100 / dt), every=round(10 / dt))
        with np.load(output / name / "fields.npz") as data:
            dt_fields[dt] = data["phi"][-1].copy()
        dt_energy[dt] = summary[name]["final_energy"]
    reference = dt_fields[0.0625]
    errors = [{"dt": dt, "rms_field_difference": float(np.sqrt(np.mean((phi-reference)**2))),
               "relative_l2_field_difference": float(np.linalg.norm(phi-reference) /
                                                      np.linalg.norm(reference)),
               "final_energy": dt_energy[dt]} for dt, phi in dt_fields.items()]
    summary["time_step_comparison"] = {"time": 100, "n": 32, "length": 32,
                                       "reference_dt": 0.0625, "comparisons": errors}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    fig, ax = plt.subplots(figsize=(5, 3.5), constrained_layout=True)
    ax.loglog([r["dt"] for r in errors[:-1]],
              [r["relative_l2_field_difference"] for r in errors[:-1]], "o-")
    ax.set(xlabel="Time step", ylabel="Relative field difference at t = 100",
           title="Reference: dt = 0.0625, fixed 32 x 32 grid")
    fig.savefig(output / "time-step-comparison.png", dpi=180)
    plt.close(fig)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("case-runs"))
    args = parser.parse_args()
    try:
        result = run_cases(args.output)
    except (ValueError, OSError, FloatingPointError) as exc:
        parser.error(str(exc))
    print(json.dumps(result["time_step_comparison"], indent=2))
