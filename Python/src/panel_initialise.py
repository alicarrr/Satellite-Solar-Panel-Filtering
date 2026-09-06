# generate the panel in solcore to get the panel type properties

import solcore
import numpy as np


def dumbDumbPanel():
    # Not using Solcore, this is just for testing. Here this is a simple panel with step functions for EQE, absorptivity, emissivity, transmissivity, and reflectivity
    wavelengths = np.arange(10, 50e3, 10)  # nm
    EQE = np.zeros_like(wavelengths)
    EQE[(wavelengths >= 400) & (wavelengths <= 800)] = 0.9
    absorptivity = np.zeros_like(wavelengths)
    absorptivity[(wavelengths >= 400) & (wavelengths <= 800)] = 0.8
    emissivity = np.zeros_like(wavelengths)
    emissivity[(wavelengths >= 400) & (wavelengths <= 800)] = 0.7
    transmissivity = np.zeros_like(wavelengths)
    transmissivity[(wavelengths >= 400) & (wavelengths <= 800)] = 0.1
    reflectivity = np.zeros_like(wavelengths)
    reflectivity[(wavelengths >= 400) & (wavelengths <= 800)] = 0.1

    panel = {
        "wavelengths": wavelengths,
        "EQE": EQE,
        "absorptivity": absorptivity,
        "emissivity": emissivity,
        "transmissivity": transmissivity,
        "reflectivity": reflectivity,
    }
    return panel