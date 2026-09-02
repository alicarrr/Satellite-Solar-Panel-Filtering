import lcm 
from include.simMessages import sim_data_t, sim_ready_t, spectrum_metadata_t, spectrum_data_t
import time
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd


receive_new_data = True  

# Write all data to a pickle file
simDataPickleFile = Path(__file__).resolve().parent / "Simulation Data" / "simulation_data.pkl"

panelLabels = ["+x", "-x", "+y", "-y", "+z", "-z"] # matches Cubesat_t struct from C++ code

sim_done = False

wavelengths = []
N_spectra = None
panelPowers = None
N_panels = None
simulation_times = []


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


def spectra_data_handler(channel, data):
    msg = spectrum_data_t.decode(data)

    global panelPowers, N_panels, N_spectra, simulation_times
    if panelPowers is None:
        N_panels = msg.N_panel
        N_spectra = msg.N_spectra
        # Array for all times
        panelPowers = []

    panelPower = np.array(msg.panelPowers)

    if panelPower.shape[0] != N_panels:
        panelPower = panelPower.T
    
    panelPowers.append(panelPower)
    simulation_times.append(msg.t)


def save_sim_data_to_pickle(pickleFile, panelPowers, simulation_times, wavelengths, panelLabels):
    panelPowersArr = np.asarray(panelPowers)
    simulation_timesArr = np.asarray(simulation_times)
    wavelengthsArr = np.asarray(wavelengths)

    if simulation_timesArr.shape[0] != panelPowersArr.shape[0]:
            raise ValueError("The number of simulation times does not match the number of power samples.")
    
    panelData = {}

    for panelIndex in range(panelPowersArr.shape[1]):
        panelData[panelLabels[panelIndex]] = pd.DataFrame(panelPowersArr[:, panelIndex, :], columns=wavelengthsArr)
        panelData[panelLabels[panelIndex]].index = simulation_timesArr
        panelData[panelLabels[panelIndex]].index.name = "simulation_time_s"
        panelData[panelLabels[panelIndex]].columns.name = "wavelength_nm"

    savedData = {"panelPowers": panelData, "simulation_times": simulation_timesArr, "wavelengths": wavelengthsArr, "panelLabels": panelLabels, "N_panels": panelPowersArr.shape[1], "N_spectra": panelPowersArr.shape[2]}
    pickleFile.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(savedData, pickleFile)
    print(f"Saved simulation data to {pickleFile}")


def load_sim_data_from_pickle():
    savedData = pd.read_pickle(simDataPickleFile)

    simulation_times = savedData["simulation_times"]
    wavelengths = np.asarray(savedData["wavelengths"])
    panelLabels = savedData["panelLabels"]
    N_panels = savedData["N_panels"]
    N_spectra = savedData["N_spectra"]
    panelPowersArr = np.stack([savedData["panelPowers"][label].to_numpy() for label in panelLabels], axis=1)

    return panelPowersArr, wavelengths, panelLabels, N_panels, N_spectra, simulation_times


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
    simulation_times_arr = np.array(simulation_times)

    # write data to pickle file
    save_sim_data_to_pickle(simDataPickleFile, panelPowers, simulation_times, wavelengths, panelLabels)
    
else:
    # Load existing data from csv files
    panelPowers_arr, wavelengths, panelLabels, N_panels, N_spectra, simulation_times = load_sim_data_from_pickle()

    # print(f'Total solar power: {panelPowers_arr.sum()}')
    panel_index = 0  # replace with the Earth-facing panel
    ir_W_m2 = panelPowers_arr[-1, panel_index, :].sum()

    print(f"Earth-facing panel IR: {ir_W_m2} W/m^2")
    print(f"Loaded simulation data from {simDataPickleFile}")

    # Plot total power for each panel over time
    for i in range(N_panels):
        plt.plot(simulation_times, panelPowers_arr[:, i, :].sum(axis=1), label=f"Panel {panelLabels[i]}")
    plt.xlabel("Simulation Time (s)")
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
    # plt.show()


    plt.figure()
    plt.plot(wavelengths, np.sum(panelPowers_arr[0, :, :], axis=0))
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Power (W)")
    # plt.xscale("log")
    plt.title("Spectrum of All Panels at Last Time Step")
    plt.legend()
    plt.grid(True, which='both', ls='--', lw=0.5)
    plt.show()