import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd


simDataPickleFile = Path(__file__).resolve().parent / "Simulation Data" / "simulation_data.pkl"


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

    plt.figure()
    for i in range(N_panels):
        plt.plot(wavelengths, panelPowersArr[-1, i, :], label=f"Panel {panelLabels[i]}")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Power (W)")
    plt.xscale("log")
    plt.title("Spectrum of Each Panel at the Last Time Sample")
    plt.legend()
    plt.grid(True, which="both", ls="--", lw=0.5)

    plt.figure()
    plt.plot(wavelengths, panelPowersArr[-1, :, :].sum(axis=0))
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Power (W)")
    plt.title("Spectrum of All Panels at the Last Time Sample")
    plt.grid(True, which="both", ls="--", lw=0.5)

    plt.show()


if __name__ == "__main__":
    plot_data()
