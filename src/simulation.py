"""
Simulation file

Note: its recommended to keep everything in one file since we copy and paste into chatgpt
a lot.
"""

import matplotlib.pyplot as plt
import time
import numpy as np
import json
import os
from dataclasses import dataclass
from typing import List

@dataclass
class TemperatureEnergyEntry:
    temperature: float
    average_energy: float

@dataclass
class TemperatureMagnetizationEntry:
    temperature: float
    average_magnetization: float

@dataclass
class SusceptibilityEntry:
    temperature: float
    susceptibility: float

@dataclass
class SpecificHeatEntry:
    temperature: float
    specific_heat: float

@dataclass
class CorrelationEntry:
    temperature: float
    distances: List[int]
    correlations: List[float]

@dataclass
class TimeCorrelationEntry:
    temperature: float
    lags: List[int]
    time_correlations: List[float]

class Simulation:
    def __init__(self, L=10, T=2.5, J=1.0, steps=1000):
        """
        Initialize the spin lattice and histories.
        L: lattice dimension
        T: temperature
        J: coupling constant
        steps: Monte Carlo steps
        """
        print("[INIT] Creating Simulation object.")
        self.L = L
        self.T = T
        self.J = J
        self.steps = steps
        # Create random spin configuration
        self.spins = np.random.choice([-1, 1], size=(L, L))
        print("[INIT] Initial spins created.")
        self.temperature_energy_history: List[TemperatureEnergyEntry] = []
        self.temperature_magnetization_history: List[TemperatureMagnetizationEntry] = []
        self.susceptibility_history: List[SusceptibilityEntry] = []
        self.specific_heat_history: List[SpecificHeatEntry] = []
        self.correlation_history: List[CorrelationEntry] = []
        self.time_correlation_history: List[TimeCorrelationEntry] = []

    def calculate_hamiltonian(self):
        """
        Calculate the total energy (Hamiltonian) of the system
        by summing neighbor interactions.
        """
        rolled_up = np.roll(self.spins, -1, axis=0)
        rolled_down = np.roll(self.spins, 1, axis=0)
        rolled_left = np.roll(self.spins, 1, axis=1)
        rolled_right = np.roll(self.spins, -1, axis=1)
        neighbors_sum = rolled_up + rolled_down + rolled_left + rolled_right
        return -self.J * np.sum(self.spins * neighbors_sum) / 2

    def calculate_magnetization(self):
        """Compute the magnetization (mean spin) of the lattice."""
        return np.mean(self.spins)

    def update_step(self):
        """
        Perform a single Metropolis update step:
        randomly pick spins, compute energy difference,
        flip spins with Boltzmann probability.
        """
        for _ in range(self.L * self.L):
            i = np.random.randint(0, self.L)
            j = np.random.randint(0, self.L)
            neighbors_sum = (
                self.spins[(i+1) % self.L, j]
                + self.spins[(i-1) % self.L, j]
                + self.spins[i, (j+1) % self.L]
                + self.spins[i, (j-1) % self.L]
            )
            dE = 2 * self.J * self.spins[i, j] * neighbors_sum
            if dE < 0 or np.random.rand() < np.exp(-dE / self.T):
                self.spins[i, j] *= -1

    def calculate_radial_correlation(self, max_distance):
        """
        Calculate the radial correlation function up to max_distance
        by averaging over all pairs of sites at distance r.
        """
        mean_spin = np.mean(self.spins)
        var_spin = np.var(self.spins)
        if var_spin == 0:
            return [0] * max_distance

        L = self.L
        corr_sums = np.zeros(max_distance)
        counts = np.zeros(max_distance)

        for i in range(L):
            for j in range(L):
                for dx in range(-max_distance, max_distance + 1):
                    for dy in range(-max_distance, max_distance + 1):
                        r = int(round(np.sqrt(dx**2 + dy**2)))
                        if 1 <= r <= max_distance:
                            neighbor_i = (i + dx) % L
                            neighbor_j = (j + dy) % L
                            corr_sums[r-1] += (self.spins[i, j] - mean_spin) * (
                                self.spins[neighbor_i, neighbor_j] - mean_spin
                            )
                            counts[r-1] += 1

        correlations = np.zeros(max_distance)
        for r in range(max_distance):
            if counts[r] > 0:
                correlations[r] = corr_sums[r] / (counts[r] * var_spin)

        return correlations.tolist()

    def calculate_time_correlation(self, snapshots, max_lag):
        """
        Compute time correlation for a set of recorded spin snapshots.
        """
        Nconfigs = len(snapshots)
        Nsites = self.L * self.L
        lags = range(max_lag + 1)
        C = np.zeros(max_lag + 1)
        # Flatten each snapshot for easy dot products
        flattened = [snap.ravel() for snap in snapshots]

        for lag in lags:
            count = 0
            corr_sum = 0.0
            for t0 in range(Nconfigs - lag):
                corr_sum += np.dot(flattened[t0], flattened[t0 + lag])
                count += 1
            # Normalization by site count and # of pairs
            C[lag] = corr_sum / (count * Nsites)

        return C.tolist()

    def simulate(
        self, 
        temperatures: List[float], 
        steps_per_temp=1000, 
        burn_in=500, 
        max_distance=None,
        compute_energy=True, 
        compute_magnetization=True, 
        compute_correlation=True,
        compute_time_correlation=False,
        max_lag=10
    ):
        """
        Iterate over given temperatures, burn in the system,
        collect data, and optionally compute correlations.
        """
        print("[SIMULATE] Starting simulation across given temperatures.")
        if max_distance is None:
            max_distance = self.L // 2
        for idx, temp in enumerate(temperatures):
            print(f"  => Setting temperature = {temp}")
            # Reinitialize spins at the first temperature in each new call
            if idx == 0:
                print("  => Reinitializing spins for new temperature range.")
                self.spins = np.random.choice([-1, 1], size=(self.L, self.L))
            self.T = temp

            print("  => Performing burn-in steps...")
            for _ in range(burn_in):
                self.update_step()

            energies = []
            magnetizations = []
            corr_accumulator = []
            snapshots = []

            print(f"  => Collecting data for {steps_per_temp} steps at T={temp}")
            for _ in range(steps_per_temp):
                self.update_step()
                if compute_energy:
                    energies.append(self.calculate_hamiltonian())
                if compute_magnetization:
                    magnetizations.append(self.calculate_magnetization())
                if compute_correlation:
                    corr_accumulator.append(self.calculate_radial_correlation(max_distance))
                if compute_time_correlation:
                    snapshots.append(self.spins.copy())

            plt.figure(figsize=(6, 6))
            plt.imshow(self.spins, cmap='gray')
            plt.title(f"Lattice Configuration at T={temp}")
            plt.colorbar(label='Spin')
            plt.savefig(f"reports/data/lattice_T{temp}.png")
            plt.close()

            # Compute average energy & specific heat
            if compute_energy:
                avg_e = np.mean(energies)
                cv = np.var(energies) / (self.T**2 * self.L**2)
                self.temperature_energy_history.append(
                    TemperatureEnergyEntry(temperature=temp, average_energy=avg_e)
                )
                self.specific_heat_history.append(
                    SpecificHeatEntry(temperature=temp, specific_heat=cv)
                )

            # Compute average magnetization & susceptibility
            if compute_magnetization:
                avg_m = np.mean(magnetizations)
                chi = np.var(magnetizations) / (self.T * self.L**2)
                self.temperature_magnetization_history.append(
                    TemperatureMagnetizationEntry(
                        temperature=temp, 
                        average_magnetization=np.abs(avg_m)
                    )
                )
                self.susceptibility_history.append(
                    SusceptibilityEntry(
                        temperature=temp, 
                        susceptibility=chi
                    )
                )

            # Average radial correlation
            if compute_correlation and corr_accumulator:
                avg_corrs = np.mean(corr_accumulator, axis=0).tolist()
                self.correlation_history.append(
                    CorrelationEntry(
                        temperature=temp,
                        distances=list(range(1, max_distance+1)),
                        correlations=avg_corrs
                    )
                )

            # Time correlation
            if compute_time_correlation and len(snapshots) > max_lag:
                time_corr = self.calculate_time_correlation(snapshots, max_lag)
                self.time_correlation_history.append(
                    TimeCorrelationEntry(
                        temperature=temp,
                        lags=list(range(max_lag + 1)),
                        time_correlations=time_corr
                    )
                )
        print("[SIMULATE] Simulation complete.")

    def _save_json(self, file_path, data, data_label):
        """
        Save data entries in JSON format.
        """
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "w") as f:
            json.dump([entry.__dict__ for entry in data], f, indent=2)
        print(f"[SAVE] {data_label} saved to {file_path}")

    def save_temperature_energy(self, temp_energy_file="reports/data/temperature_energy.json"):
        self._save_json(temp_energy_file, self.temperature_energy_history, "Temperature vs Energy data")

    def save_temperature_magnetization(self, temp_magnetization_file="reports/data/temperature_magnetization.json"):
        self._save_json(temp_magnetization_file, self.temperature_magnetization_history, "Temperature vs Magnetization data")

    def save_specific_heat(self, specific_heat_file="reports/data/specific_heat.json"):
        self._save_json(specific_heat_file, self.specific_heat_history, "Specific Heat data")

    def save_susceptibility(self, susceptibility_file="reports/data/susceptibility.json"):
        self._save_json(susceptibility_file, self.susceptibility_history, "Susceptibility data")

    def save_correlations(self, correlation_file="reports/data/correlations.json"):
        self._save_json(correlation_file, self.correlation_history, "Correlation data")

    def save_time_correlations(self, time_correlation_file="reports/data/time_correlations.json"):
        self._save_json(time_correlation_file, self.time_correlation_history, "Time Correlation data")


def plot_results(sim: Simulation, output_dir="reports/data"):
    """
    Plot and save Energy, Magnetization, Heat Capacity, and Susceptibility
    vs. Temperature.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Plot average energy vs temperature
    temps_e = [entry.temperature for entry in sim.temperature_energy_history]
    energies = [entry.average_energy for entry in sim.temperature_energy_history]
    plt.figure()
    plt.plot(temps_e, energies, 'o-')
    plt.xlabel("Temperature")
    plt.ylabel("Average Energy")
    plt.title("Energy vs. Temperature")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "energy_vs_temperature.png"))
    plt.close()

    # Plot average magnetization vs temperature
    temps_m = [entry.temperature for entry in sim.temperature_magnetization_history]
    mags = [entry.average_magnetization for entry in sim.temperature_magnetization_history]
    plt.figure()
    plt.plot(temps_m, mags, 'o-')
    plt.xlabel("Temperature")
    plt.ylabel("Average Magnetization")
    plt.title("Magnetization vs. Temperature")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "magnetization_vs_temperature.png"))
    plt.close()

    # Plot specific heat vs temperature
    temps_cv = [entry.temperature for entry in sim.specific_heat_history]
    specific_heat = [entry.specific_heat for entry in sim.specific_heat_history]
    plt.figure()
    plt.plot(temps_cv, specific_heat, 'o-')
    plt.xlabel("Temperature")
    plt.ylabel("Specific Heat")
    plt.title("Specific Heat vs. Temperature")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "specific_heat_vs_temperature.png"))
    plt.close()

    # Plot susceptibility vs temperature
    temps_chi = [entry.temperature for entry in sim.susceptibility_history]
    suscept = [entry.susceptibility for entry in sim.susceptibility_history]
    plt.figure()
    plt.plot(temps_chi, suscept, 'o-')
    plt.xlabel("Temperature")
    plt.ylabel("Susceptibility")
    plt.title("Susceptibility vs. Temperature")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, "susceptibility_vs_temperature.png"))
    plt.close()

def plot_correlations(sim: Simulation, output_dir="reports/data"):
    """
    Plot and save radial correlations for each temperature.
    """
    os.makedirs(output_dir, exist_ok=True)

    plt.figure()
    for c_entry in sim.correlation_history:
        plt.plot(
            c_entry.distances, 
            c_entry.correlations, 
            'o-', 
            label=f"T={c_entry.temperature:.2f}"
        )
    plt.xlabel("Distance")
    plt.ylabel("Correlation")
    plt.title("Radial Correlations vs. Distance")
    plt.grid(True)
    plt.legend()
    plt.savefig(os.path.join(output_dir, "radial_correlations.png"))
    plt.close()

def plot_time_correlations(sim: Simulation, output_dir="reports/data"):
    """
    Plot time correlations (autocorrelation function) as a function of lag.
    """
    os.makedirs(output_dir, exist_ok=True)
    plt.figure()
    for t_entry in sim.time_correlation_history:
        plt.plot(
            t_entry.lags, 
            t_entry.time_correlations, 
            'o-', 
            label=f"T={t_entry.temperature:.2f}"
        )
    plt.xlabel("Lag (MC Steps)")
    plt.ylabel("Time Correlation")
    plt.title("Time Correlation vs. Lag")
    plt.grid(True)
    plt.legend()
    plt.savefig(os.path.join(output_dir, "time_correlations.png"))
    plt.close()

if __name__ == "__main__":
    print("[MAIN] Running simulation")
    start_time = time.time()
    steps = 1000
    burn_in = 200

    sim = Simulation(L=10, T=1.0, J=1.0, steps=steps)

    # First: gather static data from T=0.1..5.0
    temperatures = np.linspace(0.1, 5.0, num=20)
    sim.simulate(
        temperatures,
        steps_per_temp=steps,
        burn_in=burn_in,
        compute_energy=True,
        compute_magnetization=True,
        compute_correlation=False,
        compute_time_correlation=False
    )
    sim.save_temperature_energy()
    sim.save_temperature_magnetization()
    sim.save_specific_heat()
    sim.save_susceptibility()

    # Now, collect correlation/time-correlation data at specific temps
    sim.simulate(
        [1.5, 2.5, 3.5],
        steps_per_temp=steps,
        burn_in=burn_in,
        compute_energy=False,
        compute_magnetization=False,
        compute_correlation=True,
        compute_time_correlation=True,
        max_lag=20
    )
    sim.save_correlations()
    sim.save_time_correlations()

    # Finally, create plots in the reports/data directory
    plot_results(sim, output_dir="reports/data")
    plot_correlations(sim, output_dir="reports/data")
    plot_time_correlations(sim, output_dir="reports/data")

    end_time = time.time()
    print(f"[MAIN] Execution time: {end_time - start_time:.6f} seconds")
