# All constants throughout the python scripts
import numpy as np

# Multiply degrees by this to convert to radians
DEG2RAD = np.pi / 180.0

# Multiply radians by this to convert to degrees
RAD2DEG = 180.0 / np.pi

# If there is no comment after a variable then the comments above it will appear when you hover it.

# An even better way for single constants is with constexpr, 
# it's effectively the same but instead of an unchangable variable it creates a constant value in your program.
# This means it behaves truly like a constant rather than an unchangable variable
HOUR2MINUTE = 60

m = 1.0 # m
nm = 1e-9 # m
km = 1000.0 # m
AU = 149597870700.0 # m

c = 299792458.0 # m/s

kg = 1.0 # kg
K = 1.0 # Kelvin

second = 1.0 # seconds
minute = 60.0 # seconds
hour = 60.0 * minute # seconds
day = 24.0 * hour # seconds
year = 365.24218967 * day # seconds

siderial_day = 23 * hour + 56 * minute + 4 * second # seconds
solar_day = 24 * hour # seconds

N = kg * m / (second * second) # Newton
J = N * m # Joule
W = J / second # Watt

amp = 1.0 # A
ohm = kg * m * m / (second * second * second * amp * amp) # Ohm
volt = ohm * amp # Volt
henry = ohm * second # Henry
tesla = N / (amp * m) # Tesla
farad = amp * second / volt # Farad

G = 6.67430e-11 * m*m*m / (kg * second*second) # Gravitational constant (m^3 kg^-1 s^-2)

sigma_sb = 5.670374419e-8 * W / (m * m * (K*K*K*K)) # Stefan-Boltzmann constant (W m^-2 K^-4)

M_sun = 1988500e24 * kg # Mass of the Sun (kg)
mu_sun = G * M_sun # Gravitational parameter of the Sun (m^3 s^-2)
L_sol = 3.828e26 * W # Solar luminosity (W)
R_sol = 695700e3 * m # Radius of the Sun (m)

M_earth = 5.97217e24 * kg # Mass of Earth (kg)
mu_earth = G * M_earth # Gravitational parameter of Earth (m^3 s^-2)
R_earth = 6371e3 * m # Radius of Earth (m)
albedo_earth = 0.306 # Average albedo of Earth
emissivity_IR_earth = 0.95 # Emissivity of Earth in the infrared spectrum
T_earth = 288.0 * K # Average temperature of Earth (K)
a_earth = 6378.137 * km
b_earth = 6356.7523142 * km
i_earth = 23.45 * DEG2RAD # Inclination of Earth's axis
angular_vel_earth = 2.0 * np.pi / siderial_day

lat_geomagnetic_north = 80.7 * DEG2RAD # Geomagnetic north pole latitude
lon_geomagnetic_north = -72.7 * DEG2RAD # Geomagnetic north pole longitude
magnetic_dipole_earth = 7.71e22 * amp * m * m # Earth's magnetic dipole moment (A m^2)

mu_0 = 4.0 * np.pi * 1e-7 * N / (amp * amp) # Vacuum permeability (N/A^2)

K_sun = 93.0 / W # Luminous Efficacy of the Sun (lm/W)
lux_to_footcandela = 1.0/10.7639 # Conversion factor from lux to foot-candela (lm/m^2 to lm/ft^2)

TT_delta = 69.184 * second # Difference between Terrestrial Time and Universal Time (seconds)

J2000_JD = 2451545

h = 6.62607015e-34 * J * second # Planck's constant (J s)
k = 1.380649e-23 * J / K # Boltzmann's constant (J/K)
q = 1.602176634e-19 * amp * second # Elementary charge (C)