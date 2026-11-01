# Note: this is not in the run or Cmake yet

# Pipeline
"""
--In C++
Start orbit sim with Cubesat

At each time step, get the satellite position and orientation in ECI. 

Now need to calculate the incident heat load on each panel by checking 
1. If the Sun is eclipsed by the Earth (check if the line from the satellite to the Sun intersects the Earth sphere)
2. If the panel is facing the sun (dot product of panel normal and sun vector)
3. If not eclipsed, which parts of Earth are visible to the panel, and what is the average albedo of those parts of Earth (use HAMSTER data to get the albedo spectrum at the visible lat/lon points)
4. If in view of Earth, the fraction visible (of the total Earth surface area) will emit Earth IR radiation, which can be approximated as a greybody

Output should be a time series of the spectral incident heat load on each panel (send to LCM channel for use in Python)
----

--In Python
Read the spectral heat load time series from LCM for each panel. Save LCM Message to a file.

Panel model (from package) including:
- spectral absorptivity, emissivity, transmissivity, and reflectivity
- Calculation with SolCore of the Quantum Efficiency spectrum of the panel (regardless of incident heat load)
- Potential for Computation of efficiency and other performance metrics
- These performance metrics can be temperature dependent (like total power output, panel efficiency, etc.)

Temperature calculation (to feed into the panel model) - likely just basic Stefan Boltzmann law with root finding (as it's in a vacuum, the only heat loss is through radiation)

Filter model (Casadi symbolic) to calculate the effect of a filter on the panel performance in each spectral property (absorptivity, emissivity, transmissivity, and reflectivity).

Non-linear optimiser (Casadi)
- Define a cost function. 
    - For a panel type take in the quantum efficiency, spectral properties, and the incident heat load time series.
    - Initialise the optimisation for a filter spectrum (which can be a symbolic vector of the same length as the wavelength range of the incident heat load time series)
    - Initially, just maximise for power output without temperature dependence.
    - Can then penalise thermal load (to reduce panel temperature). This would require calculating temperature of the panel within the optimisation, so the relationship between temperature and properties would need to be defined analytically or through a lookup table.
- Find the optimal filter spectra to maximize panel performance (efficiency, power output, etc.)


"""