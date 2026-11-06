import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

# latex formatting
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


simDataPickleFile = Path(__file__).resolve().parent / "Simulation Data" / "simulation_data.pkl"
spectrum_snapshots = [(20.0, "sunlit"), (65.0, "eclipse")]  # minutes from simulation start


def load_sim_data_from_pickle(pickleFile):
    savedData = pd.read_pickle(pickleFile)

    simulation_times = np.asarray(savedData["simulation_times"])
    wavelengths = np.asarray(savedData["wavelengths"])
    print(f'Max wavelength: {wavelengths.max()} nm, Min wavelength: {wavelengths.min()} nm')
    panelLabels = savedData["panelLabels"]
    N_panels = savedData["N_panels"]
    N_spectra = savedData["N_spectra"]
    panelPowersArr = np.stack([savedData["panelPowers"][label].to_numpy() for label in panelLabels], axis=1)
    panelTypes = savedData["panelTypes"]

    return panelPowersArr, wavelengths, panelLabels, panelTypes, N_panels, N_spectra, simulation_times


def plot_data():
    panelPowersArr, wavelengths, panelLabels, panelTypes, N_panels, N_spectra, simulation_times = load_sim_data_from_pickle(simDataPickleFile)
    print(f"Panel Types: {panelTypes}")

    print(f"Loaded simulation data from {simDataPickleFile}")
    print(f"Data shape: {panelPowersArr.shape} = [time, panel, wavelength]")

    plt.figure()
    for i in range(N_panels):
        plt.plot(simulation_times, panelPowersArr[:, i, :].sum(axis=1), label=f"Panel {panelLabels[i]}")
    plt.xlabel("Simulation Time (s)")
    plt.ylabel("Total Power (W)")
    plt.title("Power Incident on Each Panel Over Time")
    plt.legend()
    plt.grid(True, ls="--", lw=0.5)

    # Plot figures of the spectral load at a given time, one in eclipse and one in sunlight (based on estimates of time)
    figures = simDataPickleFile.parent / "figures"

    figures.mkdir(exist_ok=True)

    for requested_minutes, state in spectrum_snapshots:
        requested_time = simulation_times[0] + 60.0 * requested_minutes
        if not simulation_times[0] <= requested_time <= simulation_times[-1]:
            raise ValueError(f"The {state} snapshot at {requested_minutes:g} min is outside the saved simulation. Use the full-orbit simulation_data.pkl.")
        sample = int(np.argmin(np.abs(simulation_times - requested_time)))
        actual_minutes = (simulation_times[sample] - simulation_times[0]) / 60.0
        fig, axes = plt.subplots(3, 2, figsize=(6.3, 6.0), sharex=True)
        for i, (ax, label) in enumerate(zip(axes.flat, panelLabels)):
            ax.plot(wavelengths, panelPowersArr[sample, i, :], linewidth=0.7)
            ax.set(xscale="log", xlim=(wavelengths.min(), wavelengths.max()))
            ax.set_ylim(bottom=0)
            ax.set_title(f"({chr(97 + i)}) ${label}$", loc="left")
            if i >= 4:
                ax.set_xlabel(r"Wavelength, $\lambda$ (nm)")
            ax.grid(alpha=0.25)
        fig.supylabel(r"Incident power per bin, $\dot Q_{\mathrm{in},i}(\lambda_k)$ (W)")

        # fig.suptitle(f"{state.capitalize()} snapshot, $t={actual_minutes:.2f}$ min")
        fig.tight_layout()
        for extension in ("pdf", "png"):
            fig.savefig(figures / f"orbital_incident_spectrum_{state}_{actual_minutes:.2f}_min.{extension}", dpi=200, bbox_inches="tight")
        print(f"Saved {state} spectra at t = {actual_minutes:.2f} min (sample {sample}).")

    plt.figure()
    plt.plot(wavelengths, panelPowersArr[-1, :, :].sum(axis=0))
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Power (W)")
    plt.title("Spectrum of All Panels at the Last Time Sample")
    plt.grid(True, which="both", ls="--", lw=0.5)

    plt.show()


if __name__ == "__main__":
    plot_data()
