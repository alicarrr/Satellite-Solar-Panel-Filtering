import lcm 
from include.simMessages import sim_data_t, sim_ready_t, spectrum_metadata_t, spectrum_data_t
import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd

####################
# NOTE: Need to change data save to use time instead of time step index
####################
receive_new_data = True  
# Write all data to a pickle file
simDataPickleFile = Path(__file__).resolve().parent / "Simulation Data" / "simulation_data.pkl"

panelLabels = ["+x", "-x", "+y", "-y", "+z", "-z"] # matches Cubesat_t struct from C++ code

global sim_done
sim_done = False

global wavelengths, N_spectra
wavelengths = []

global panelPowers, N_panels
panelPowers = None
N_panels = None


def sim_finished_handler(channel, data):
    msg = sim_ready_t.decode(data)
    if (msg.is_ready == True):
        print("PY: Received Simulation Finished")
        global sim_done
        sim_done = True

def sim_data_handler(channel, data):
    msg = sim_data_t.decode(data)
    print(f"Received data on channel: {channel}")
    print(msg.t)
    print(f"pos = {msg.pos_sat_eci}")

def spectra_metadata_handler(channel, data):
    msg = spectrum_metadata_t.decode(data)

    print("PY: Received Metadata")
    global wavelengths, N_spectra
    wavelengths = np.array(msg.wavelengths)
    N_spectra = wavelengths.shape[0]
    # print(wavelengths.shape)


def spectra_data_handler(channel, data):
    msg = spectrum_data_t.decode(data)

    global panelPowers, N_panels
    if (panelPowers is None):
        N_panels = msg.N_panel
        N_spectra = msg.N_spectra

        #One time step of panel powers
        panelPower = np.array(msg.panelPowers)
        if (panelPower.shape[0] != N_panels):
            panelPower = panelPower.T

        # Array for all times
        panelPowers = []
        panelPowers.append(panelPower)

    else:
        panelPower = np.array(msg.panelPowers)
        if (panelPower.shape[0] != N_panels):
            panelPower = panelPower.T

        # Array for all times
        panelPowers.append(panelPower)

def save_sim_data_to_pickle(pickleFile, panelPowers, wavelengths, panelLabels):
    panelPowersArr = np.asarray(panelPowers)
    wavelengthsArr = np.asarray(wavelengths)
    panelData = {}

    for panelIndex in range(panelPowersArr.shape[1]):
        panelData[panelLabels[panelIndex]] = pd.DataFrame(panelPowersArr[:, panelIndex, :], columns=wavelengthsArr)
        panelData[panelLabels[panelIndex]].index.name = "time_step"
        panelData[panelLabels[panelIndex]].columns.name = "wavelength_nm"

    savedData = {"panelPowers": panelData, "wavelengths": wavelengthsArr, "panelLabels": panelLabels, "N_panels": panelPowersArr.shape[1], "N_spectra": panelPowersArr.shape[2]}
    pickleFile.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(savedData, pickleFile)
    print(f"Saved simulation data to {pickleFile}")


def load_sim_data_from_pickle():
    savedData = pd.read_pickle(simDataPickleFile)
    wavelengths = np.asarray(savedData["wavelengths"])
    panelLabels = savedData["panelLabels"]
    N_panels = savedData["N_panels"]
    N_spectra = savedData["N_spectra"]
    panelPowersArr = np.stack([savedData["panelPowers"][label].to_numpy() for label in panelLabels], axis=1)
    return panelPowersArr, wavelengths, panelLabels, N_panels, N_spectra


time.sleep(2)

if receive_new_data:
    lc = lcm.LCM()

    msg = sim_ready_t()
    msg.is_ready = True

    lc.publish("ListenerReady", msg.encode())



    metadata_sub = lc.subscribe("SpectraMetaData", spectra_metadata_handler)
    spectra_sub = lc.subscribe("SpectraData", spectra_data_handler)
    sim_done_sub = lc.subscribe("SimStatus", sim_finished_handler)

    try:
        while (not sim_done):
            # keep checking for new messages
            lc.handle()
    except KeyboardInterrupt:
        pass

    # If you make it here then simulation is finished


    # Indexed by [time][panel num][wavelength index]
    panelPowers_arr = np.array(panelPowers)

    # write data to pickle file
    save_sim_data_to_pickle(simDataPickleFile, panelPowers, wavelengths, panelLabels)
    
else:
    # Load existing data from csv files
    panelPowers_arr, wavelengths, panelLabels, N_panels, N_spectra = load_sim_data_from_pickle()

    print(f'Total solar power: {panelPowers_arr.sum()}')

    print(f"Loaded simulation data from {simDataPickleFile}")

    # Plot total power for each panel over time
    for i in range(N_panels):
        plt.plot(panelPowers_arr[:, i, :].sum(axis=1), label=f"Panel {panelLabels[i]}")
    plt.xlabel("Time")
    plt.ylabel("Total Power (W)")
    plt.title("Power Output of Solar Panels Over Time")
    plt.legend()

    # Plot the spectrum for each panel at the last time step
    plt.figure()
    for i in range(N_panels):
        plt.plot(wavelengths, panelPowers_arr[0, i, :], label=f"Panel {panelLabels[i]}")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Power (W)")
    plt.xscale("log")
    plt.title("Spectrum of Solar Panels at Last Time Step")
    plt.legend()
    plt.grid(True, which='both', ls='--', lw=0.5)
    plt.show()