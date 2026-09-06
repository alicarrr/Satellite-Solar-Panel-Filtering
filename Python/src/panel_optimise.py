import casadi as ca
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path


PLANCK_CONSTANT = 6.62607015e-34
SPEED_OF_LIGHT = 299792458.0
ELEMENTARY_CHARGE = 1.602176634e-19


def load_simulation_data(pickleFile, manualPanelTypes=None):
    savedData = pd.read_pickle(pickleFile)

    times = np.asarray(savedData["simulation_times"], dtype=float)
    wavelengths = np.asarray(savedData["wavelengths"], dtype=float)
    panelLabels = savedData["panelLabels"]
    panelPowers = np.stack([savedData["panelPowers"][label].to_numpy() for label in panelLabels], axis=1)

    if manualPanelTypes is not None:
        panelTypes = np.asarray(manualPanelTypes, dtype=int)
    elif "panelTypes" in savedData:
        panelTypes = np.asarray(savedData["panelTypes"], dtype=int)
    else:
        raise ValueError("The pickle file does not contain panelTypes, so manualPanelTypes must be supplied.")

    if panelPowers.shape[1] != panelTypes.size:
        raise ValueError("The number of panel types does not match the number of panels.")
    if panelPowers.shape[0] != times.size or panelPowers.shape[2] != wavelengths.size:
        raise ValueError("The saved panel-power dimensions do not match the time and wavelength arrays.")
    if times.size < 2 or np.any(np.diff(times) <= 0.0):
        raise ValueError("Simulation times must contain at least two strictly increasing values.")
    if wavelengths.size < 2 or np.any(np.diff(wavelengths) <= 0.0):
        raise ValueError("Wavelengths must contain at least two strictly increasing values.")

    return panelPowers, times, wavelengths, panelLabels, panelTypes


def integrate_incident_energy(panelPowers, times):
    N_panels = panelPowers.shape[1]
    N_wavelengths = panelPowers.shape[2]
    incidentEnergy = np.zeros((N_panels, N_wavelengths))

    for timeIndex in range(times.size - 1):
        dt = times[timeIndex + 1] - times[timeIndex]
        averagePower = 0.5*(panelPowers[timeIndex] + panelPowers[timeIndex + 1])
        incidentEnergy += averagePower*dt

    return incidentEnergy


def evaluate_spectrum(value, wavelengths, name):
    if callable(value):
        spectrum = np.asarray(value(wavelengths), dtype=float)
    elif np.isscalar(value):
        spectrum = np.full(wavelengths.shape, float(value))
    else:
        spectrum = np.asarray(value, dtype=float)

    if spectrum.shape != wavelengths.shape:
        raise ValueError(f"{name} must return one value for every wavelength.")
    if np.any(~np.isfinite(spectrum)):
        raise ValueError(f"{name} contains a non-finite value.")

    return spectrum


def build_substrate_data(wavelengths, substrateDefinitions):
    substrateData = {}
    wavelengths_m = wavelengths*1.0e-9

    for panelType, definition in substrateDefinitions.items():
        reflectivity = evaluate_spectrum(definition["reflectivity"], wavelengths, f"type {panelType} reflectivity")
        transmissivity = evaluate_spectrum(definition["transmissivity"], wavelengths, f"type {panelType} transmissivity")
        absorptivity = evaluate_spectrum(definition["absorptivity"], wavelengths, f"type {panelType} absorptivity")

        if np.any(reflectivity < 0.0) or np.any(transmissivity < 0.0) or np.any(absorptivity < 0.0):
            raise ValueError(f"Type {panelType} has a negative optical property.")
        if np.any(reflectivity > 1.0) or np.any(transmissivity > 1.0) or np.any(absorptivity > 1.0):
            raise ValueError(f"Type {panelType} has an optical property greater than one.")
        if not np.allclose(reflectivity + transmissivity + absorptivity, 1.0, atol=1.0e-6):
            raise ValueError(f"Type {panelType} reflectivity, transmissivity and absorptivity must sum to one.")

        if "quantumEfficiency" in definition:
            quantumEfficiency = evaluate_spectrum(definition["quantumEfficiency"], wavelengths, f"type {panelType} quantum efficiency")
            operatingVoltage = float(definition["operatingVoltage"])

            if operatingVoltage <= 0.0:
                raise ValueError(f"Type {panelType} operating voltage must be positive.")

            electricalEfficiency = quantumEfficiency*ELEMENTARY_CHARGE*operatingVoltage*wavelengths_m/(PLANCK_CONSTANT*SPEED_OF_LIGHT)
        else:
            quantumEfficiency = np.zeros_like(wavelengths)
            electricalEfficiency = np.zeros_like(wavelengths)

        if np.any(quantumEfficiency < 0.0) or np.any(quantumEfficiency > 1.0):
            raise ValueError(f"Type {panelType} quantum efficiency must be between zero and one.")
        if np.any(electricalEfficiency > absorptivity + 1.0e-10):
            raise ValueError(f"Type {panelType} electrical efficiency exceeds its optical absorptivity.")

        heatAbsorptivity = np.maximum(absorptivity - electricalEfficiency, 0.0)
        substrateData[int(panelType)] = {"reflectivity": reflectivity, "transmissivity": transmissivity, "absorptivity": absorptivity, "quantumEfficiency": quantumEfficiency, "electricalEfficiency": electricalEfficiency, "heatAbsorptivity": heatAbsorptivity}

    return substrateData


def interpolate_knot_values(wavelengths, knotWavelengths, knotValues):
    interpolatedValues = []

    for wavelength in wavelengths:
        if wavelength <= knotWavelengths[0]:
            interpolatedValues.append(knotValues[0])
        elif wavelength >= knotWavelengths[-1]:
            interpolatedValues.append(knotValues[-1])
        else:
            rightIndex = int(np.searchsorted(knotWavelengths, wavelength))
            leftIndex = rightIndex - 1
            fraction = (wavelength - knotWavelengths[leftIndex])/(knotWavelengths[rightIndex] - knotWavelengths[leftIndex])
            value = (1.0 - fraction)*knotValues[leftIndex] + fraction*knotValues[rightIndex]
            interpolatedValues.append(value)

    return ca.vertcat(*interpolatedValues)


def create_knot_wavelengths(wavelengths, maximumStep_nm=10.0):
    if maximumStep_nm <= 0.0:
        raise ValueError("The maximum coating-knot step must be positive.")

    start = wavelengths[0]
    stop = wavelengths[-1]
    N_fullSteps = int(np.floor((stop - start)/maximumStep_nm))
    knotWavelengths = start + maximumStep_nm*np.arange(N_fullSteps + 1)

    if stop - knotWavelengths[-1] > 1.0e-9:
        knotWavelengths = np.append(knotWavelengths, stop)
    else:
        knotWavelengths[-1] = stop

    return knotWavelengths


def calculate_performance(incidentEnergy, wavelengths, panelTypes, substrateData, coatingGroupByType, coatingSpectra):
    groupNames = list(dict.fromkeys(coatingGroupByType[int(panelType)] for panelType in panelTypes))
    groupResults = {}

    for groupName in groupNames:
        groupResults[groupName] = {"incidentEnergySpectrum_J": np.zeros(wavelengths.size), "electricalEnergySpectrum_J": np.zeros(wavelengths.size), "heatEnergySpectrum_J": np.zeros(wavelengths.size)}

    totalElectricalEnergy = 0.0
    totalHeatEnergy = 0.0

    for panelIndex in range(panelTypes.size):
        panelType = int(panelTypes[panelIndex])
        groupName = coatingGroupByType[panelType]
        substrate = substrateData[panelType]
        coating = coatingSpectra[groupName]

        for wavelengthIndex in range(wavelengths.size):
            incident = incidentEnergy[panelIndex, wavelengthIndex]
            coatingTransmissivity = coating["transmissivity"][wavelengthIndex]
            coatingAbsorptivity = coating["absorptivity"][wavelengthIndex]
            electricalEfficiency = substrate["electricalEfficiency"][wavelengthIndex]
            heatAbsorptivity = substrate["heatAbsorptivity"][wavelengthIndex]

            electrical = incident*coatingTransmissivity*electricalEfficiency
            heat = incident*(coatingAbsorptivity + coatingTransmissivity*heatAbsorptivity)

            totalElectricalEnergy += electrical
            totalHeatEnergy += heat
            groupResults[groupName]["incidentEnergySpectrum_J"][wavelengthIndex] += incident
            groupResults[groupName]["electricalEnergySpectrum_J"][wavelengthIndex] += electrical
            groupResults[groupName]["heatEnergySpectrum_J"][wavelengthIndex] += heat

    return {"electricalEnergy_J": totalElectricalEnergy, "heatEnergy_J": totalHeatEnergy, "groups": groupResults}


def optimise_coatings(incidentEnergy, wavelengths, panelTypes, substrateData, coatingGroupByType, knotWavelengths, minimumElectricalFraction=0.999, smoothnessWeight=1.0e-4):
    knotWavelengths = np.asarray(knotWavelengths, dtype=float)

    if np.any(incidentEnergy < 0.0):
        raise ValueError("Incident energy cannot be negative.")
    if minimumElectricalFraction < 0.0 or minimumElectricalFraction > 1.0:
        raise ValueError("minimumElectricalFraction must be between zero and one.")
    if knotWavelengths.size < 3 or np.any(np.diff(knotWavelengths) <= 0.0):
        raise ValueError("At least three strictly increasing coating knots are required.")
    if knotWavelengths[0] > wavelengths[0] or knotWavelengths[-1] < wavelengths[-1]:
        raise ValueError("The coating knots must cover the full wavelength range.")

    for panelType in np.unique(panelTypes):
        if int(panelType) not in coatingGroupByType:
            raise ValueError(f"No coating group was supplied for panel type {panelType}.")

    opti = ca.Opti()
    coatings = {}
    smoothness = 0
    groupNames = list(dict.fromkeys(coatingGroupByType[int(panelType)] for panelType in panelTypes))
    baselineCoatings = {}

    for groupName in groupNames:
        baselineCoatings[groupName] = {"reflectivity": np.zeros(wavelengths.size), "transmissivity": np.ones(wavelengths.size), "absorptivity": np.zeros(wavelengths.size)}

    baselinePerformance = calculate_performance(incidentEnergy, wavelengths, panelTypes, substrateData, coatingGroupByType, baselineCoatings)

    for groupName in groupNames:
        reflectivityKnots = opti.variable(knotWavelengths.size)
        transmissivityKnots = opti.variable(knotWavelengths.size)

        opti.subject_to(opti.bounded(0.0, reflectivityKnots, 1.0))
        opti.subject_to(opti.bounded(0.0, transmissivityKnots, 1.0))
        opti.subject_to(reflectivityKnots + transmissivityKnots <= 1.0)

        for knotIndex in range(1, knotWavelengths.size - 1):
            reflectivityCurvature = reflectivityKnots[knotIndex + 1] - 2.0*reflectivityKnots[knotIndex] + reflectivityKnots[knotIndex - 1]
            transmissivityCurvature = transmissivityKnots[knotIndex + 1] - 2.0*transmissivityKnots[knotIndex] + transmissivityKnots[knotIndex - 1]
            smoothness += reflectivityCurvature**2 + transmissivityCurvature**2

        opti.set_initial(reflectivityKnots, 0.1)
        opti.set_initial(transmissivityKnots, 0.8)

        reflectivity = interpolate_knot_values(wavelengths, knotWavelengths, reflectivityKnots)
        transmissivity = interpolate_knot_values(wavelengths, knotWavelengths, transmissivityKnots)
        absorptivity = 1.0 - reflectivity - transmissivity
        coatings[groupName] = {"reflectivity": reflectivity, "transmissivity": transmissivity, "absorptivity": absorptivity, "reflectivityKnots": reflectivityKnots, "transmissivityKnots": transmissivityKnots}

    electricalEnergy = 0
    absorbedHeat = 0

    for panelIndex in range(panelTypes.size):
        panelType = int(panelTypes[panelIndex])

        if panelType not in substrateData:
            raise ValueError(f"No substrate data was supplied for panel type {panelType}.")
        if panelType not in coatingGroupByType:
            raise ValueError(f"No coating group was supplied for panel type {panelType}.")

        substrate = substrateData[panelType]
        coating = coatings[coatingGroupByType[panelType]]

        for wavelengthIndex in range(wavelengths.size):
            incident = float(incidentEnergy[panelIndex, wavelengthIndex])
            coatingTransmissivity = coating["transmissivity"][wavelengthIndex]
            coatingAbsorptivity = coating["absorptivity"][wavelengthIndex]
            electricalEfficiency = float(substrate["electricalEfficiency"][wavelengthIndex])
            heatAbsorptivity = float(substrate["heatAbsorptivity"][wavelengthIndex])

            electricalEnergy += incident*coatingTransmissivity*electricalEfficiency
            absorbedHeat += incident*(coatingAbsorptivity + coatingTransmissivity*heatAbsorptivity)

    totalIncidentEnergy = float(np.sum(incidentEnergy))

    if totalIncidentEnergy <= 0.0:
        raise ValueError("The total incident energy must be positive.")

    maximumElectricalEnergy = baselinePerformance["electricalEnergy_J"]

    if maximumElectricalEnergy > 0.0:
        opti.subject_to(electricalEnergy >= minimumElectricalFraction*maximumElectricalEnergy)

    objective = absorbedHeat/totalIncidentEnergy + smoothnessWeight*smoothness

    opti.minimize(objective)
    opti.solver("ipopt", {"expand": True}, {"max_iter": 2000, "tol": 1.0e-8})
    solution = opti.solve()

    coatingResults = {}

    for groupName, coating in coatings.items():
        coatingResults[groupName] = {"reflectivity": np.asarray(solution.value(coating["reflectivity"])).squeeze(), "transmissivity": np.asarray(solution.value(coating["transmissivity"])).squeeze(), "absorptivity": np.asarray(solution.value(coating["absorptivity"])).squeeze()}

    optimisedPerformance = calculate_performance(incidentEnergy, wavelengths, panelTypes, substrateData, coatingGroupByType, coatingResults)
    results = {"coatings": coatingResults, "baseline": baselinePerformance, "optimised": optimisedPerformance, "objective": float(solution.value(objective)), "minimumElectricalFraction": minimumElectricalFraction}
    return results


def example_solar_quantum_efficiency(wavelengths):
    quantumEfficiency = np.zeros_like(wavelengths)
    quantumEfficiency[(wavelengths >= 400.0) & (wavelengths <= 1100.0)] = 0.85
    return quantumEfficiency


def constant_spectrum(value):
    return lambda wavelengths: np.full(wavelengths.shape, value, dtype=float)


def plot_coating_results(wavelengths, panelTypes, substrateData, coatingGroupByType, results):
    for groupName, coating in results["coatings"].items():
        baseline = results["baseline"]["groups"][groupName]
        optimised = results["optimised"]["groups"][groupName]
        groupPanelTypes = [int(panelType) for panelType in np.unique(panelTypes) if coatingGroupByType[int(panelType)] == groupName]
        figure, axes = plt.subplots(6, 1, figsize=(10, 15), sharex=True)

        axes[0].plot(wavelengths, baseline["incidentEnergySpectrum_J"])
        axes[0].set_ylabel("Incident Energy\nper Bin (J)")
        axes[0].set_title(f"Inputs and Result for the {groupName} Coating")

        for panelType in groupPanelTypes:
            substrate = substrateData[panelType]
            quantumEfficiency = substrate["quantumEfficiency"]
            axes[1].fill_between(wavelengths, 0.0, quantumEfficiency, alpha=0.25, label=f"Type {panelType} quantum efficiency")
            axes[1].plot(wavelengths, quantumEfficiency, lw=1.0)
            axes[1].plot(wavelengths, substrate["electricalEfficiency"], label=f"Type {panelType} electrical efficiency")
            axes[2].plot(wavelengths, substrate["absorptivity"], label=f"Type {panelType} absorptivity")
            axes[2].plot(wavelengths, substrate["reflectivity"], label=f"Type {panelType} reflectivity")
            axes[2].plot(wavelengths, substrate["transmissivity"], label=f"Type {panelType} transmissivity")
            axes[2].plot(wavelengths, substrate["heatAbsorptivity"], ls="--", label=f"Type {panelType} non-electrical absorption")

        axes[1].set_ylabel("QE or Electrical\nEfficiency")
        axes[1].set_ylim(0.0, 1.0)
        axes[1].legend()

        axes[2].set_ylabel("Underlying Panel\nProperty")
        axes[2].set_ylim(0.0, 1.0)
        axes[2].legend()

        axes[3].plot(wavelengths, coating["transmissivity"], label="Transmissivity")
        axes[3].plot(wavelengths, coating["reflectivity"], label="Reflectivity")
        axes[3].plot(wavelengths, coating["absorptivity"], label="Absorptivity")
        axes[3].set_ylabel("Coating Property")
        axes[3].set_ylim(0.0, 1.0)
        axes[3].legend()

        axes[4].plot(wavelengths, baseline["electricalEnergySpectrum_J"], ls="--", label="No coating")
        axes[4].plot(wavelengths, optimised["electricalEnergySpectrum_J"], label="Optimised coating")
        axes[4].set_ylabel("Electrical Energy\nper Bin (J)")
        axes[4].legend()

        axes[5].plot(wavelengths, baseline["heatEnergySpectrum_J"], ls="--", label="No coating")
        axes[5].plot(wavelengths, optimised["heatEnergySpectrum_J"], label="Optimised coating")
        axes[5].set_xlabel("Wavelength (nm)")
        axes[5].set_ylabel("Heat Energy\nper Bin (J)")
        axes[5].legend()

        for axis in axes:
            axis.set_xscale("log")
            axis.grid(True, which="both", ls="--", lw=0.5)

        figure.tight_layout()

    figure, axes = plt.subplots(1, 2, figsize=(9, 4))
    labels = ["No coating", "Optimised coating"]
    electricalValues = [results["baseline"]["electricalEnergy_J"], results["optimised"]["electricalEnergy_J"]]
    heatValues = [results["baseline"]["heatEnergy_J"], results["optimised"]["heatEnergy_J"]]

    axes[0].bar(labels, electricalValues)
    axes[0].set_ylabel("Electrical Energy (J)")
    axes[0].set_title("Electrical-Energy Comparison")
    axes[0].grid(True, axis="y", ls="--", lw=0.5)

    axes[1].bar(labels, heatValues)
    axes[1].set_ylabel("Absorbed Heat (J)")
    axes[1].set_title("Heat-Energy Comparison")
    axes[1].grid(True, axis="y", ls="--", lw=0.5)

    figure.tight_layout()

    plt.show()


def main():
    pickleFile = Path(__file__).resolve().parent / "Simulation Data" / "simulation_data.pkl"

    # Replace these example types with the actual type of each panel. Set this to None once panelTypes is saved in the pickle file.
    manualPanelTypes = np.array([0, 0, 0, 0, 1, 1])

    panelPowers, times, wavelengths, panelLabels, panelTypes = load_simulation_data(pickleFile, manualPanelTypes)
    incidentEnergy = integrate_incident_energy(panelPowers, times)

    # These are example underlying-surface properties and must be replaced with measured or selected spectra.
    substrateDefinitions = {
        0: {"reflectivity": constant_spectrum(0.10), "transmissivity": constant_spectrum(0.00), "absorptivity": constant_spectrum(0.90), "quantumEfficiency": example_solar_quantum_efficiency, "operatingVoltage": 0.50},
        1: {"reflectivity": constant_spectrum(0.85), "transmissivity": constant_spectrum(0.00), "absorptivity": constant_spectrum(0.15)}
    }

    substrateData = build_substrate_data(wavelengths, substrateDefinitions)

    # All panels of type 0 share one coating, and all panels of type 1 share another coating.
    coatingGroupByType = {0: "solar", 1: "aluminium"}
    knotWavelengths = create_knot_wavelengths(wavelengths, maximumStep_nm=10.0)

    print(f"Using {knotWavelengths.size} coating knots with a maximum step of {np.max(np.diff(knotWavelengths))} nm")

    results = optimise_coatings(incidentEnergy, wavelengths, panelTypes, substrateData, coatingGroupByType, knotWavelengths, minimumElectricalFraction=0.999, smoothnessWeight=1.0e-4)

    print(f"No-coating electrical energy: {results['baseline']['electricalEnergy_J']} J")
    print(f"Optimised electrical energy: {results['optimised']['electricalEnergy_J']} J")
    print(f"No-coating absorbed heat: {results['baseline']['heatEnergy_J']} J")
    print(f"Optimised absorbed heat: {results['optimised']['heatEnergy_J']} J")
    print(f"Minimum retained electrical fraction: {results['minimumElectricalFraction']}")
    print(f"Objective: {results['objective']}")

    plot_coating_results(wavelengths, panelTypes, substrateData, coatingGroupByType, results)


if __name__ == "__main__":
    main()
