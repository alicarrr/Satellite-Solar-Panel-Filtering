import lcm 
from include.simMessages import sim_data_t, sim_ready_t, spectrum_metadata_t, spectrum_data_t
import time
import numpy as np
import matplotlib.pyplot as plt

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

    



time.sleep(2)

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

# Plot total power for each panel over time
for i in range(N_panels):
    plt.plot(panelPowers_arr[:, i, :].sum(axis=1), label=f"Panel {i+1}")

plt.xlabel("Time")
plt.ylabel("Total Power")
plt.title("Power Output of Solar Panels Over Time")
plt.legend()
plt.show()

print(panelPowers_arr.shape)