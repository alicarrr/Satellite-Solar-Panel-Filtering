import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib as mpl

import pvlib as pv

# Plotting formatting with LaTex and auto adding grid lines.
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "axes.grid": True,
    "grid.alpha": 0.5,
    "grid.linestyle": "--",
    "grid.linewidth": 0.5,
})



wavelength_arr = np.linspace(0, 1100, int(1e4))


# Power output of sun emissivity function and panel absortivity function, as well as scale factor (dot product of sun vector and panel normal vector), include area?

def sun_emissivity_function(wavelength):
    # Test with square wave emissivity function.
    return np.where(wavelength < 500, 0.9, 0.5)


def panel_absorptivity_function(wavelength):        
    # Always between 0 and 1, test with square wave function.
    return np.where(wavelength < 700, 0.8, 0.3)


def total_power_output(wavelength_arr):
     # For now, not using scale factor, assumes sun vector and panel normal vector are perfectly aligned,
    # and area is 1 m^2.
    def integrand(wavelength):
        return sun_emissivity_function(wavelength) * panel_absorptivity_function(wavelength)
    
    total_power = integrate.quad(integrand, float(wavelength_arr[0]), float(wavelength_arr[-1]))[0]
    
    return total_power

print(f'{total_power_output(wavelength_arr)} W/m^2')

# Sample plots
# Absorption spectrum of the panel
plt.figure()
plt.plot(wavelength_arr, panel_absorptivity_function(wavelength_arr))
plt.title("Absorption Spectrum of the Panel")
plt.xlabel(r"$\lambda$ (nm)")
plt.ylabel("Absorptivity")
plt.xlim(wavelength_arr[0], wavelength_arr[-1])
plt.savefig("absorption_spectrum.png", dpi=300)



# Emission spectrum of the sun
plt.figure()
plt.plot(wavelength_arr, sun_emissivity_function(wavelength_arr))
plt.title("Emission Spectrum of the Sun")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Emissivity")
plt.xlim(wavelength_arr[0], wavelength_arr[-1])
plt.savefig("emission_spectrum.png", dpi=300)


# Power output over each wavelength, which is the product of the above two functions.
plt.figure()
plt.plot(wavelength_arr, sun_emissivity_function(wavelength_arr) * panel_absorptivity_function(wavelength_arr))
plt.title("Power Output at Each Wavelength")
plt.xlabel("Wavelength (nm)")
plt.ylabel("Power Output")
plt.xlim(wavelength_arr[0], wavelength_arr[-1])
plt.text(600, 0.5, f"Total Power Output (W/m$^2$): {total_power_output(wavelength_arr):.2f}", fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
plt.savefig("power_output.png", dpi=300)

# plt.show()
