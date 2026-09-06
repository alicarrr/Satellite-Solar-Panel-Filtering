import lcm
from include.simMessages import sim_ready_t, spectrum_metadata_t, spectrum_data_t
import time
import numpy as np
from pathlib import Path
import pandas as pd


simDataPickleFile = Path(__file__).resolve().parent / "Simulation Data" / "simulation_data.pkl"
panelLabels = ["+x", "-x", "+y", "-y", "+z", "-z"] # matches Cubesat_t struct from C++ code

sim_done = False
wavelengths = []
N_spectra = None
panelPowers = []
N_panels = None
simulation_times = []


def sim_finished_handler(channel, data):
    msg = sim_ready_t.decode(data)
    if msg.is_ready:
        print("PY: Received Simulation Finished")
        global sim_done
        sim_done = True


def spectra_metadata_handler(channel, data):
    msg = spectrum_metadata_t.decode(data)
    print("PY: Received Metadata")

    global wavelengths, N_spectra
    wavelengths = np.asarray(msg.wavelengths)
    N_spectra = wavelengths.shape[0]


def spectra_data_handler(channel, data):
    msg = spectrum_data_t.decode(data)

    global N_panels, N_spectra
    if N_panels is None:
        N_panels = msg.N_panel
        N_spectra = msg.N_spectra

    panelPower = np.asarray(msg.panelPowers)
    if panelPower.shape == (N_spectra, N_panels):
        panelPower = panelPower.T
    elif panelPower.shape != (N_panels, N_spectra):
        panelPower = panelPower.reshape(N_panels, N_spectra)

    panelPowers.append(panelPower)
    simulation_times.append(msg.t)


def save_sim_data_to_pickle(pickleFile, panelPowers, simulation_times, wavelengths, panelLabels):
    panelPowersArr = np.asarray(panelPowers)
    simulation_timesArr = np.asarray(simulation_times)
    wavelengthsArr = np.asarray(wavelengths)

    if panelPowersArr.ndim != 3:
        raise ValueError("Panel power data must have dimensions [time, panel, wavelength].")
    if simulation_timesArr.shape[0] != panelPowersArr.shape[0]:
        raise ValueError("The number of simulation times does not match the number of power samples.")
    if wavelengthsArr.shape[0] != panelPowersArr.shape[2]:
        raise ValueError("The number of wavelengths does not match the spectrum size.")
    if len(panelLabels) != panelPowersArr.shape[1]:
        raise ValueError("The number of panel labels does not match the number of panels.")

    panelData = {}
    for panelIndex in range(panelPowersArr.shape[1]):
        panelData[panelLabels[panelIndex]] = pd.DataFrame(panelPowersArr[:, panelIndex, :], columns=wavelengthsArr, index=simulation_timesArr)
        panelData[panelLabels[panelIndex]].index.name = "simulation_time_s"
        panelData[panelLabels[panelIndex]].columns.name = "wavelength_nm"

    savedData = {"panelPowers": panelData, "simulation_times": simulation_timesArr, "wavelengths": wavelengthsArr, "panelLabels": panelLabels, "N_panels": panelPowersArr.shape[1], "N_spectra": panelPowersArr.shape[2]}
    pickleFile.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(savedData, pickleFile)
    print(f"Saved simulation data to {pickleFile}")


def main():
    lc = lcm.LCM()
    lc.subscribe("SpectraMetaData", spectra_metadata_handler)
    lc.subscribe("SpectraData", spectra_data_handler)
    lc.subscribe("SimStatus", sim_finished_handler)

    time.sleep(2)

    msg = sim_ready_t()
    msg.is_ready = True
    lc.publish("ListenerReady", msg.encode())

    try:
        while not sim_done:
            lc.handle()
    except KeyboardInterrupt:
        print("PY: Listener stopped by user")

    if not panelPowers:
        print("No spectrum data was received, so no pickle file was written.")
        return

    save_sim_data_to_pickle(simDataPickleFile, panelPowers, simulation_times, wavelengths, panelLabels)


if __name__ == "__main__":
    main()
