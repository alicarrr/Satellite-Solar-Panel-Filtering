"""CubeSat sensitivity study. Place beside optimise_filters.py and run normally.

Edit the settings below. The default study has thirteen cases, each with five
starts: one nominal guess and four random guesses. The original optimiser's
equations and data files are reused, including its electrical approximation.

Each run creates a new folder with case pickles, terminal logs, summary CSVs
and PDF/PNG figures. Each pickle includes all starts and the best successful
solution. Case settings are also saved as JSON files. Set n_workers to run
several cases at once. Each case's five starts run sequentially, reusing its
built optimiser. The console shows progress and full solver output goes directly
to each case's text file while it runs.

Set reoptimise = False to regenerate tables and plots from saved pickles.
That mode reads the saved results without needing the original model files.
Plot formatting uses LaTeX and Computer Modern, so a working LaTeX installation
is needed when creating figures.
"""

import csv
import itertools
import json
import os
import pickle
import subprocess
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Latex format in plots
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 12,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12
})


# =============================================================================
# 1. Settings to edit
# =============================================================================
here = Path(__file__).resolve().parent
output_root = here / "Optimisation Results" / "CubeSat Sensitivity Study"

reoptimise = False
# None selects the newest study with saved pickles when reoptimise is False.
# To select another run, set this to its directory (absolute or relative to here).
saved_study_directory = None

n_starts = 3                    # total starts: one nominal + four random
n_workers = 2                   # cases running at once, 1 runs them sequentially
threads_per_worker = 1          # numerical-library threads in each process
random_seed = 0                # same random starts across cases for comparison
verbose_ipopt = True           # full iteration output, saved as it prints
emission_averaging_temperature_k = 273.15  # fixed numerical reference for all cases

baseline = {
    "internal_heat_load_w": 10.0,       # total CubeSat internal heat [W]
    "t_set": 273.15,                   # temperature target for all faces [K]
    "w_elec": 100.0,                   # w_E, electrical energy weight
    "w_temp": 1.0,                     # normalised temperature penalty weight
    "w_smooth": 1.0e-3,                # w_S, filter roughness weight
    "incident_power_scale": 1.0       # multiplies solar + albedo + Earth IR
}

# Include the baseline value in each list. A singleton means no sweep.
# Defaults give 13 unique cases / 65 solves in one_at_a_time mode.
# Vary each weight relative to the others. Scaling all three equally leaves
# the mathematical optimum unchanged. Compare energy, tracking and roughness.
# w_smooth = 0 removes the roughness penalty, the slope constraints still apply.
parameter_values = {
    "internal_heat_load_w": [0.0, 10.0, 20.0],
    "t_set": [253.15, 273.15, 293.15, 313.15],
    "w_elec": [10.0, 100.0, 1000.0],
    "w_temp": [0.1, 1.0, 10.0],
    "w_smooth": [0.0, 1.0e-3, 1.0, 1000.0],  # wide range: baseline is a small tie-breaker
    "incident_power_scale": [1.0]       # e.g. [0.75, 1.0, 1.25]
}
sweep_mode = "one_at_a_time"            # or "full_factorial" for all combinations

make_case_plots = True
make_comparison_plots = True
figure_formats = ("pdf", "png")

parameter_labels = {
    "internal_heat_load_w": "Total internal heat load [W]",
    "t_set": "Target temperature [K]",
    "w_elec": r"Electrical weight $w_E$",
    "w_temp": r"Temperature weight $w_T$",
    "w_smooth": r"Smoothness weight $w_S$",
    "incident_power_scale": "Incident radiation multiplier"
}


# =============================================================================
# 2. Case generation and file helpers
# =============================================================================
def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2), encoding="utf-8")


def save_pickle(path, data):
    # Replace only after the write completes, so interrupted writes are obvious.
    path = Path(path)
    temporary = path.with_suffix(".pkl.tmp")
    with temporary.open("wb") as stream:
        pickle.dump(data, stream, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(path)


def load_pickle(path):
    with Path(path).open("rb") as stream:
        return pickle.load(stream)


def make_cases():
    if set(parameter_values) != set(baseline):
        raise ValueError("baseline and parameter_values must have the same keys.")
    if set(baseline) != set(parameter_labels):
        raise ValueError("Supported sweep parameters are listed in parameter_labels.")
    if sweep_mode not in ("one_at_a_time", "full_factorial"):
        raise ValueError("sweep_mode must be one_at_a_time or full_factorial.")
    if not isinstance(n_starts, int) or n_starts < 1:
        raise ValueError("n_starts must be a positive integer.")
    if not isinstance(n_workers, int) or n_workers < 1:
        raise ValueError("n_workers must be a positive integer.")
    if not isinstance(threads_per_worker, int) or threads_per_worker < 1:
        raise ValueError("threads_per_worker must be a positive integer.")
    if not np.isfinite(emission_averaging_temperature_k) or emission_averaging_temperature_k <= 0:
        raise ValueError("emission_averaging_temperature_k must be finite and positive.")
    for name, values in parameter_values.items():
        if not values or baseline[name] not in values:
            raise ValueError(f"Include the baseline value in parameter_values[{name!r}].")
        for value in values:
            if not np.isfinite(value) or value < 0:
                raise ValueError(f"{name} values must be finite and non-negative.")
            if name in ("t_set", "incident_power_scale") and value == 0:
                raise ValueError(f"{name} must be strictly positive.")

    settings = [dict(baseline)]
    if sweep_mode == "one_at_a_time":
        for name, values in parameter_values.items():
            settings.extend({**baseline, name: float(value)} for value in values)
    else:
        names = list(parameter_values)
        for values in itertools.product(*parameter_values.values()):
            settings.append(dict(zip(names, values)))

    cases = []
    for parameters in settings:
        parameters = {k: float(v) for k, v in parameters.items()}
        if not any(parameters[name] > 0 for name in ("w_elec", "w_temp", "w_smooth")):
            raise ValueError("Each case must have at least one positive objective weight.")
        if any(case["parameters"] == parameters for case in cases):
            continue
        changed = [k for k in baseline if parameters[k] != baseline[k]]
        if not changed:
            label = "baseline"
        elif len(changed) == 1:
            name = changed[0]
            label = f"{name}_{parameters[name]:g}"
        else:
            label = "combined"
        cases.append({"case_id": f"{len(cases):03d}_{label}", "parameters": parameters, "changed_parameters": changed})
    return cases


def write_csv(path, rows):
    if not rows:
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


# =============================================================================
# 3. Physical metrics, including plots-only mode
# =============================================================================
def calculate_metrics(payload):
    # Energies use the optimiser's time steps, with zero weight on the last sample.
    # Whole-satellite efficiency includes radiation on all six faces. PV efficiency
    # includes only the four PV faces. Both use incident power before the filter.
    # Tracking RMSE weights all panel/time samples equally, matching the objective.
    labels = payload["panel_labels"]
    t = np.asarray(payload["t"])
    dt = np.r_[np.diff(t), 0.0]
    duration = float(t[-1] - t[0])
    pv = np.asarray(payload["material"]) == payload["study"]["solar_material"]
    settings = payload["settings"]
    incident = np.asarray(payload["incident_energy_per_panel_J"])

    def measure(power, temperature):
        energy = float(np.sum(power * dt[None, :]))
        error = temperature - settings["t_set"]
        return {
            "electrical_energy_J": energy,
            "mean_electrical_power_W": energy / duration,
            "peak_electrical_power_W": float(np.max(power.sum(axis=0))),
            "whole_satellite_efficiency_pct": 100.0 * energy / incident.sum(),
            "pv_efficiency_pct": (100.0 * np.sum(power[pv] * dt[None, :]) / incident[pv].sum() if incident[pv].sum() > 0 else np.nan),
            "tracking_rmse_K": float(np.sqrt(np.mean(error ** 2))),
            "tracking_mae_K": float(np.mean(np.abs(error))),
            "tracking_max_abs_error_K": float(np.max(np.abs(error))),
            "pv_tracking_rmse_K": float(np.sqrt(np.mean(error[pv] ** 2))),
            "aluminium_tracking_rmse_K": float(np.sqrt(np.mean(error[~pv] ** 2))),
            "minimum_temperature_K": float(temperature.min()),
            "maximum_temperature_K": float(temperature.max())
        }

    filtered_t = np.stack([payload["T"][lab] for lab in labels])
    metrics = measure(np.asarray(payload["electrical_power_W"]), filtered_t)
    reference = measure(np.asarray(payload["unfiltered_electrical_power_W"]), np.asarray(payload["unfiltered_temperature_K"]))
    metrics.update({"unfiltered_" + k: v for k, v in reference.items()})
    gain = metrics["electrical_energy_J"] - reference["electrical_energy_J"]
    metrics["electrical_energy_gain_J"] = gain
    metrics["electrical_energy_gain_pct"] = (100 * gain / reference["electrical_energy_J"] if reference["electrical_energy_J"] > 0 else np.nan)

    # Exactly the same roughness normalisation as build_problem(), Section 5.6.
    roughness = sum(np.sum(np.diff(f[name]) ** 2) for f in payload["filters"].values() for name in ("alpha", "rho", "tau"))
    roughness /= len(payload["filters"]) * len(payload["lam_nm"])
    metrics["filter_roughness"] = float(roughness)
    metrics["electrical_cost"] = (-settings["w_elec"] * metrics["electrical_energy_J"] / incident.sum())
    metrics["temperature_cost"] = (settings["w_temp"] * (metrics["tracking_rmse_K"] / settings["t_ref_scale"]) ** 2)
    metrics["smoothness_cost"] = settings["w_smooth"] * roughness
    metrics["objective"] = float(payload["J"])
    metrics["objective_reconstructed"] = sum(metrics[k] for k in ("electrical_cost", "temperature_cost", "smoothness_cost"))
    return metrics


# =============================================================================
# 4. Run one case using the original optimiser
# =============================================================================
def run_case(config_path):
    config_path = Path(config_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    output_path = config_path.with_suffix(".pkl")
    started = time.perf_counter()
    result = {"study": config, "success": False, "start_results": []}
    try:
        # Import here so plots-only mode does not need CasADi or the model files.
        import optimise_filters as optimiser

        print(f"\nCASE {config['case_id']}: {config['parameters']}", flush=True)

        parameters = config["parameters"]
        optimiser.internal_heat_load_w = parameters["internal_heat_load_w"]
        optimiser.w_elec = parameters["w_elec"]
        optimiser.w_temp = parameters["w_temp"]
        optimiser.w_smooth = parameters["w_smooth"]

        # Pass grid settings explicitly: the original function's default
        # arguments were bound at import time.
        data = optimiser.load_sim_data(optimiser.simulation_pickle_file, optimiser.n_filter_bands, optimiser.band_spacing)
        expected = {"+x": optimiser.solar, "-x": optimiser.solar, "+y": optimiser.solar, "-y": optimiser.solar, "+z": optimiser.aluminium, "-z": optimiser.aluminium}
        actual = dict(zip(data.labels, map(int, data.material)))
        if data.n_p != 6 or actual != expected:
            raise ValueError(f"Expected four side PV panels and two aluminium ends, got {actual}.")

        # Scale both the working-grid power and native absorption weights.
        # This keeps the original spectrum/time history and needs no large copy.
        scale = parameters["incident_power_scale"]
        data.q_in *= scale
        data.native_q_spectrum *= scale
        # In the original code t_set also weights emissivity band averages.
        # Hold that numerical approximation fixed when changing the target.
        optimiser.t_set = config["emission_averaging_temperature_K"]
        props = optimiser.build_panel_properties(data)
        optimiser.t_set = parameters["t_set"]
        options = dict(optimiser.ipopt_options)
        if not config["verbose_ipopt"]:
            options.update({"ipopt.print_level": 0, "ipopt.sb": "yes"})
        prob = optimiser.build_problem(data, props, solver_options=options)

        config["solar_material"] = int(optimiser.solar)
        config["solver_options"] = options

        rng = np.random.default_rng(config["random_seed"])
        solutions = []
        for start in range(config["n_starts"]):
            print(f"\n--- Start {start + 1}/{config['n_starts']} (index {start}) ---", flush=True)
            start_time = time.perf_counter()
            record = {"start": start, "success": False}
            try:
                guess = optimiser.make_initial_guess(prob, None if start == 0 else rng)
                solution = optimiser.solve_once(prob, guess, start=start)
                # vars() gives ordinary dictionaries, no optimiser class is saved.
                record.update(vars(solution))
                solution.success = bool(solution.success and np.isfinite(solution.J))
                record["success"] = solution.success
                solutions.append(solution)
                stats = prob.opti.stats()
                record["return_status"] = stats.get("return_status", "unknown")
                record["iterations"] = stats.get("iter_count")
            except Exception:
                record["message"] = traceback.format_exc()
                print(record["message"], flush=True)
            record["elapsed_s"] = time.perf_counter() - start_time
            result["start_results"].append(record)
            print(f"Start {start}: success={record['success']}, J={record.get('J', float('nan')):.8g}, wall time={record['elapsed_s']:.1f} s", flush=True)
            # Checkpoint after every start, even if all attempts so far failed.
            result["elapsed_s"] = time.perf_counter() - started
            save_pickle(output_path, result)

        successful = [s for s in solutions if s.success]
        if not successful:
            raise RuntimeError("None of the starts succeeded. All attempts remain in this case's pickle.")
        best = min(successful, key=lambda s: s.J)
        result.update(optimiser.solution_payload(best, solutions, prob, props))
        result.update({"success": True, "best_start": best.start, "n_successful_starts": len(successful), "incident_power_per_panel_W": data.q_in.sum(axis=2), "incident_energy_per_panel_J": np.sum(data.q_in * data.dt[None, :, None], axis=(1, 2))})
        result["metrics"] = calculate_metrics(result)
        optimiser.report(best, prob)
    except Exception:
        result["success"] = False
        result["error"] = traceback.format_exc()
        print(result["error"], flush=True)
    result["elapsed_s"] = time.perf_counter() - started
    save_pickle(output_path, result)
    print(f"Saved {output_path}", flush=True)
    return 0 if result["success"] else 1


# =============================================================================
# 5. Run the study and save terminal output
# =============================================================================
def optimise_study():
    cases = make_cases()
    if not (here / "optimise_filters.py").is_file():
        raise FileNotFoundError("Place this script beside optimise_filters.py.")
    run_directory = output_root / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    run_directory.mkdir(parents=True, exist_ok=False)
    manifest = {"baseline": baseline, "parameter_values": parameter_values, "sweep_mode": sweep_mode, "cases": cases, "n_starts": n_starts, "n_workers": n_workers, "threads_per_worker": threads_per_worker, "random_seed": random_seed, "emission_averaging_temperature_K": emission_averaging_temperature_k}
    write_json(run_directory / "manifest.json", manifest)

    # Set these before starting Python so the child imports NumPy with this limit.
    # One thread per process lets separate cases use the available CPU cores.
    environment = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        environment[name] = str(threads_per_worker)

    running = {}
    next_case = 0
    with (run_directory / "terminal_output.txt").open("w", encoding="utf-8") as combined:
        def echo(line):
            print(line, end="", flush=True)
            combined.write(line)
            combined.flush()

        echo(f"Study: {run_directory}\n{len(cases)} cases, {n_starts} starts per case, up to {n_workers} cases at once.\n")
        try:
            while next_case < len(cases) or running:
                # Start another case whenever a worker becomes available.
                while next_case < len(cases) and len(running) < n_workers:
                    case = cases[next_case]
                    next_case += 1
                    config = {**case, "n_starts": n_starts, "random_seed": random_seed, "verbose_ipopt": verbose_ipopt, "threads_per_worker": threads_per_worker, "emission_averaging_temperature_K": emission_averaging_temperature_k}
                    config_path = run_directory / (case["case_id"] + ".json")
                    write_json(config_path, config)
                    command = [sys.executable, "-u", str(Path(__file__).resolve()), "--run-case", str(config_path)]
                    # A separate file captures Python and native IPOPT output.
                    with config_path.with_suffix(".txt").open("w", encoding="utf-8") as log:
                        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=environment)
                    running[case["case_id"]] = (process, config_path)
                    echo(f"Started {case['case_id']} (PID {process.pid}). Log: {config_path.with_suffix('.txt').name}\n")

                for case_id, (process, config_path) in list(running.items()):
                    return_code = process.poll()
                    if return_code is None:
                        continue
                    status = "completed" if return_code == 0 else f"FAILED (exit {return_code})"
                    echo(f"Case {case_id}: {status}.\n")
                    # Also collect the full logs in one file, grouped by case.
                    combined.write(f"\n--- Full log: {case_id} ---\n")
                    with config_path.with_suffix(".txt").open(encoding="utf-8", errors="replace") as log:
                        for line in log:
                            combined.write(line)
                    combined.flush()
                    del running[case_id]
                if running:
                    time.sleep(0.2)
        finally:
            # Stop active cases if the user interrupts the study.
            for process, config_path in running.values():
                if process.poll() is None:
                    process.terminate()
            for process, config_path in running.values():
                process.wait()
    return run_directory


# =============================================================================
# 6. Tables and figures from saved pickles
# =============================================================================
def read_study(directory):
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    results, summary, starts = [], [], []
    for case in manifest["cases"]:
        path = directory / (case["case_id"] + ".pkl")
        payload = load_pickle(path) if path.exists() else {"success": False, "error": "No saved result."}
        row = {"case_id": case["case_id"], **case["parameters"], "success": payload["success"], "elapsed_s": payload.get("elapsed_s", np.nan), "n_starts_attempted": len(payload.get("start_results", []))}
        if payload["success"]:
            payload["metrics"] = calculate_metrics(payload)
            row.update(payload["metrics"])
            row.update({"best_start": payload["best_start"], "n_successful_starts": payload["n_successful_starts"]})
            results.append(payload)
        else:
            row["error"] = payload.get("error", "Run interrupted before a successful result was saved.")
        summary.append(row)
        for record in payload.get("start_results", []):
            start_row = {"case_id": case["case_id"]}
            for name in ("start", "success", "J", "e_elec", "efficiency", "t_track", "return_status", "iterations", "elapsed_s", "message"):
                start_row[name] = record.get(name)
            starts.append(start_row)
    write_csv(directory / "summary.csv", summary)
    write_csv(directory / "starts.csv", starts)
    lines = ["Case | energy [J] | mean power [W] | PV efficiency [%] | tracking RMSE [K] | J"]
    for row in summary:
        if row["success"]:
            lines.append(f"{row['case_id']} | {row['electrical_energy_J']:.4f} | {row['mean_electrical_power_W']:.4f} | {row['pv_efficiency_pct']:.4f} | {row['tracking_rmse_K']:.4f} | {row['objective']:.6f}")
        else:
            lines.append(f"{row['case_id']} | FAILED (see case log/pickle)")
    text = "\n".join(lines) + "\n"
    print(text)
    (directory / "summary.txt").write_text(text, encoding="utf-8")
    return manifest, results


def make_figures(directory, manifest, results):
    # Plot all figures
    figures = directory / "figures"
    figures.mkdir(exist_ok=True)

    def save(fig, name):
        fig.tight_layout()
        for extension in figure_formats:
            fig.savefig(figures / f"{name}.{extension}", dpi=200, bbox_inches="tight")
        plt.close(fig)

    if make_case_plots:
        for payload in results:
            tag = payload["study"]["case_id"]
            times = np.asarray(payload["t"])
            minutes = (times - times[0]) / 60.0

            nominal = payload["study"]["parameters"] == manifest["baseline"]
            material = np.asarray(payload["material"])
            solar_material = payload["study"]["solar_material"]
            filtered_colour, unfiltered_colour = "#0072B2", "#D55E00"
            filtered_temperature = np.asarray([payload["T"][label] for label in payload["panel_labels"]])
            unfiltered_temperature = np.asarray(payload["unfiltered_temperature_K"])

            temperature_limits = {}

            for group in np.unique(material):
                values = np.r_[filtered_temperature[material == group].ravel(), unfiltered_temperature[material == group].ravel(), payload["settings"]["t_set"]]
                padding = max(1.0, 0.05 * np.ptp(values))
                temperature_limits[group] = (values.min() - padding, values.max() + padding)
            fig, axes = plt.subplots(3, 2, figsize=(6.5, 5.3), sharex=True)

            for i, (ax, label) in enumerate(zip(axes.flat, payload["panel_labels"])):
                ax.plot(minutes, unfiltered_temperature[i], "--", color=unfiltered_colour, label="Unfiltered")
                ax.plot(minutes, filtered_temperature[i], color=filtered_colour, label="Filtered")
                ax.axhline(payload["settings"]["t_set"], color="0.25", ls=":", label="Target")
                if "t_bounds" in payload and unfiltered_temperature[i].max() > payload["t_bounds"][1]:
                    ax.axhline(payload["t_bounds"][1], color="0.55", ls="-.", lw=0.7, label="Upper bound")
                    # ax.axhline(payload["t_bounds"][0], color="0.55", ls="-.", lw=0.7)
                name = "(silicon PV)" if material[i] == solar_material else "(aluminium)"
                ax.set_title(f"({chr(97 + i)}) ${label}$ {name}", loc="left")
                ax.set(xlim=(minutes[0], minutes[-1]), ylim=temperature_limits[material[i]])

                if i % 2 == 0:
                    ax.set_ylabel("Temperature (K)")
                if i >= 4:
                    ax.set_xlabel("Time (min)")
                ax.grid(alpha=0.25)
            legend = {}

            for ax in axes.flat:
                handles, labels = ax.get_legend_handles_labels()
                legend.update(zip(labels, handles))
            fig.legend(legend.values(), legend.keys(), loc="lower center", bbox_to_anchor=(0.5, 0.97), ncol=4, frameon=False)
            save(fig, tag + "_temperature")

            fig, axes = plt.subplots(2, 1, figsize=(6.5, 4.7), sharex=True)

            for key, colour, style, label in (("unfiltered_electrical_power_W", unfiltered_colour, "--", "Unfiltered"), ("electrical_power_W", filtered_colour, "-", "Filtered")):
                power = np.asarray(payload[key]).sum(axis=0)
                # Use the optimiser's left-endpoint interval weights. The final
                # sample has no following interval and contributes no energy.
                energy_kj = np.r_[0.0, np.cumsum(power[:-1] * np.diff(times))] / 1000.0
                axes[0].plot(minutes, power, style, color=colour, label=label)
                axes[1].plot(minutes, energy_kj, style, color=colour, label=label)
            axes[0].set_title("(a) Total electrical power", loc="left")
            axes[0].set_ylabel("Electrical power (W)")
            axes[0].legend(ncol=2)
            axes[1].set_title("(b) Cumulative electrical energy", loc="left")
            axes[1].set(xlabel="Time (min)", ylabel="Electrical energy (kJ)")
            for ax in axes:
                ax.set_xlim(minutes[0], minutes[-1])
                ax.set_ylim(bottom=0)
                ax.grid(alpha=0.25)
            save(fig, tag + "_power")

            wavelength_nm = np.asarray(payload["lam_nm"])
            wavelength_edges = np.asarray(payload.get("lam_edges_nm", payload["lam_nm"]))
            fig, axes = plt.subplots(2, 1, figsize=(6.5, 4.7), sharex=True)
            for ax, (filter_material, curves) in zip(axes, sorted(payload["filters"].items())):
                for name, symbol, style in (("alpha", r"$\alpha_f$", "-"), ("rho", r"$\rho_f$", "--"), ("tau", r"$\tau_f$", "-.")):
                    ax.plot(wavelength_nm, curves[name], style, label=symbol)
                ax.set(xscale="log", xlim=(wavelength_edges[0], wavelength_edges[-1]), ylim=(0, 1), ylabel="Optical fraction", title="PV filter" if filter_material == solar_material else "Aluminium filter")
                ax.grid(alpha=0.25)
                ax.legend(ncol=3, loc="center left")
            axes[-1].set_xlabel(r"Wavelength, $\lambda$ (nm)")
            save(fig, tag + "_filters")

            if nominal:
                coefficients = (("alpha", r"$\alpha_f$", "-"), ("rho", r"$\rho_f$", "--"), ("tau", r"$\tau_f$", "-."))
                for filter_material, curves in sorted(payload["filters"].items()):
                    fig, ax = plt.subplots(figsize=(6.5, 2.7))
                    for name, symbol, style in coefficients:
                        ax.plot(wavelength_nm, curves[name], style, label=symbol)
                    ax.set(xscale="log", xlim=(wavelength_edges[0], wavelength_edges[-1]), ylim=(0, 1), xlabel=r"Wavelength, $\lambda$ (nm)", ylabel="Optical fraction")
                    ax.grid(alpha=0.25)
                    ax.legend(ncol=3, loc="center left")
                    suffix = "_pv_filter" if filter_material == solar_material else "_aluminium_filter"
                    save(fig, tag + suffix)

                power_fig, power_ax = plt.subplots(figsize=(6.5, 2.7))
                energy_fig, energy_ax = plt.subplots(figsize=(6.5, 2.7))
                power_styles = (("unfiltered_electrical_power_W", unfiltered_colour, "--", "Unfiltered"), ("electrical_power_W", filtered_colour, "-", "Filtered"))
                for key, colour, style, label in power_styles:
                    power = np.asarray(payload[key]).sum(axis=0)
                    energy_kj = np.r_[0.0, np.cumsum(power[:-1] * np.diff(times))] / 1000.0
                    power_ax.plot(minutes, power, style, color=colour, label=label)
                    energy_ax.plot(minutes, energy_kj, style, color=colour, label=label)

                power_panels = ((power_fig, power_ax, "_electrical_power", "Electrical power (W)"), (energy_fig, energy_ax, "_cumulative_energy", "Electrical energy (kJ)"))
                for fig, ax, suffix, ylabel in power_panels:
                    ax.set(xlim=(minutes[0], minutes[-1]), xlabel="Time (min)", ylabel=ylabel)
                    ax.set_ylim(bottom=0)
                    ax.grid(alpha=0.25)
                    ax.legend(ncol=2)
                    save(fig, tag + suffix)

            # The saved per-panel totals are common to the nominal comparison.
            # Export them once, without requiring the original simulation file.
            if nominal and "incident_power_per_panel_W" in payload:
                incident = np.asarray(payload["incident_power_per_panel_W"])
                solar_peak = incident[material == solar_material].max() if np.any(material == solar_material) else 0.0
                fig, axes = plt.subplots(3, 2, figsize=(6.3, 6.0), sharex=True)
                for i, (ax, label) in enumerate(zip(axes.flat, payload["panel_labels"])):
                    ax.plot(minutes, incident[i], color=filtered_colour, lw=0.75)
                    peak = solar_peak if material[i] == solar_material else incident[i].max()
                    ax.set(xlim=(minutes[0], minutes[-1]), ylim=(0, max(1.05 * peak, 1.0e-12)))
                    name = "(silicon PV)" if material[i] == solar_material else "(aluminium)"
                    ax.set_title(f"({chr(97 + i)}) ${label}$ {name}", loc="left")
                    if i >= 4:
                        ax.set_xlabel("Time (min)")
                    ax.grid(alpha=0.25)
                fig.supylabel("Incident Power (W)")
                save(fig, "orbital_incident_power")

            # fig, axes = plt.subplots(3, 2, figsize=(9, 8), sharex=True, sharey=True)
            # for i, (ax, label) in enumerate(zip(axes.flat, payload["panel_labels"])):
            #     ax.plot(minutes, payload["T"][label], label="Filtered")
            #     ax.plot(minutes, payload["unfiltered_temperature_K"][i], "--", label="Unfiltered")
            #     ax.axhline(payload["settings"]["t_set"], color="k", ls=":", label="Target")
            #     ax.set(title=label, xlabel="Time [min]", ylabel="Temperature [K]")
            #     ax.grid(alpha=0.25)
            # axes.flat[0].legend(fontsize=8)
            # save(fig, tag + "_temperature")

            # fig, ax = plt.subplots(figsize=(6.5, 4))
            # ax.plot(minutes, np.asarray(payload["electrical_power_W"]).sum(axis=0), label="Filtered")
            # ax.plot(minutes, np.asarray(payload["unfiltered_electrical_power_W"]).sum(axis=0), "--", label="Unfiltered")
            # ax.set(xlabel="Time [min]", ylabel="Total electrical power [W]")
            # ax.grid(alpha=0.25)
            # ax.legend()
            # save(fig, tag + "_power")

            # fig, axes = plt.subplots(2, 1, figsize=(6.5, 6), sharex=True)
            # for ax, (material, curves) in zip(axes, sorted(payload["filters"].items())):
            #     for name, symbol, style in (("alpha", r"$\alpha_f$", "-"), ("rho", r"$\rho_f$", "--"), ("tau", r"$\tau_f$", "-.")):
            #         ax.plot(payload["lam_nm"], curves[name], style, label=symbol)
            #     ax.set(xscale="log", ylim=(0, 1), ylabel="Optical fraction", title="PV filter" if material == payload["study"]["solar_material"] else "Aluminium filter")
            #     ax.grid(alpha=0.25)
            #     ax.legend(ncol=3)
            # axes[-1].set_xlabel("Wavelength [nm]")
            # save(fig, tag + "_filters")

    if make_comparison_plots:
        # Hold other parameters at baseline values for these comparison plots.
        # Combined cases are still included in the tables and individual plots.
        baseline = manifest["baseline"]
        plot_metrics = (("electrical_energy_J", "Electrical energy [J]"), ("pv_efficiency_pct", r"PV conversion efficiency [\%]"), ("tracking_rmse_K", "Tracking RMSE [K]"), ("objective", "Case objective J"))
        for name in manifest["parameter_values"]:
            selected = []
            for payload in results:
                parameters = payload["study"]["parameters"]
                other_changes = [key for key in baseline if key != name and parameters[key] != baseline[key]]
                if not other_changes:
                    selected.append(payload)
            selected.sort(key=lambda p: p["study"]["parameters"][name])
            if len(selected) < 2:
                continue
            x = [p["study"]["parameters"][name] for p in selected]
            # Include the unweighted roughness when testing cost-function weights.
            case_metrics = plot_metrics
            if name in ("w_elec", "w_temp", "w_smooth"):
                case_metrics += (("filter_roughness", "Filter roughness"), ("mean_electrical_power_W", "Mean electrical power [W]"))
            n_columns = len(case_metrics) // 2
            fig, axes = plt.subplots(2, n_columns, figsize=(4.5 * n_columns, 6))
            for ax, (metric, label) in zip(axes.flat, case_metrics):
                ax.plot(x, [p["metrics"][metric] for p in selected], "o-", label="Filtered")
                reference_key = "unfiltered_" + metric
                if reference_key in selected[0]["metrics"]:
                    ax.plot(x, [p["metrics"][reference_key] for p in selected], "s--", label="Unfiltered")
                ax.axvline(baseline[name], color="0.6", lw=1, ls=":")
                ax.set(xlabel=parameter_labels[name], ylabel=label)
                if name in ("w_elec", "w_temp", "w_smooth"):
                    positive_weights = [value for value in x if value > 0]
                    if min(x) > 0:
                        ax.set_xscale("log")
                    elif positive_weights:
                        ax.set_xscale("symlog", linthresh=min(positive_weights))
                ax.grid(alpha=0.25)
            axes.flat[0].legend()
            save(fig, "sensitivity_" + name)


def main():
    # Private worker entry point, normal use is just `python this_file.py`.
    if len(sys.argv) > 1:
        if len(sys.argv) == 3 and sys.argv[1] == "--run-case":
            return run_case(sys.argv[2])
        raise ValueError("Edit Section 1 and run this script without arguments.")
    if reoptimise:
        directory = optimise_study()
    elif saved_study_directory is not None:
        directory = Path(saved_study_directory)
        if not directory.is_absolute():
            directory = here / directory
    else:
        saved = sorted(p.parent for p in output_root.glob("*/manifest.json") if any(p.parent.glob("*.pkl")))
        if not saved:
            raise FileNotFoundError(f"No saved study in {output_root}. Set saved_study_directory explicitly.")
        directory = saved[-1]
    manifest, results = read_study(directory)
    if make_case_plots or make_comparison_plots:
        make_figures(directory, manifest, results)
    print(f"Results: {directory}")
    if not results:
        print("No successful cases. See summary.csv and the individual logs.")
        return 1
    if len(results) != len(manifest["cases"]):
        print(f"{len(manifest['cases']) - len(results)} cases failed or are incomplete. See summary.csv.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
