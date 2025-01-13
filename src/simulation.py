import faulthandler
faulthandler.enable()

import numpy as np
import pandas as pd
import os
from numba import njit
import json
import typer

@njit
def precompute_neighbors_numba(L):
    indices = np.arange(L)
    neighbors_up = np.roll(indices, -1)
    neighbors_down = np.roll(indices, 1)
    return neighbors_up, neighbors_down

@njit
def update_step_numba(L, spins, neighbors_up, neighbors_down, J, T, h_values, step):
    """
    Perform a single update step using Metropolis dynamics.
    """
    for _ in range(L * L):
        i = np.random.randint(0, L)  # Random index within bounds
        j = np.random.randint(0, L)

        # Compute neighbors sum
        neighbors_sum = (
            spins[neighbors_up[i], j]
            + spins[neighbors_down[i], j]
            + spins[i, neighbors_up[j]]
            + spins[i, neighbors_down[j]]
        )

        # External field
        field = h_values[step, i, j] if h_values.shape[0] > step else 0

        # Energy difference if we flip this spin
        dE = 2.0 * J * spins[i, j] * (neighbors_sum + field)

        # Metropolis acceptance rule
        if dE <= 0.0:
            spins[i, j] *= -1
        else:
            if np.random.rand() < np.exp(-dE / T):
                spins[i, j] *= -1

@njit
def run_simulation(L, steps, burn_in_steps, save_every, spins, neighbors_up, neighbors_down, J, T, h_values):
    """
    Run the simulation with a burn-in phase followed by data collection.
    """
    # ---------------------------
    # Burn-in phase (no data saved)
    # ---------------------------
    for _ in range(burn_in_steps):
        update_step_numba(L, spins, neighbors_up, neighbors_down, J, T, h_values, 0)
    
    # ---------------------------
    # Main simulation (with data)
    # ---------------------------
    data = []
    for step in range(steps):
        # Use current 'step' for indexing h_values if needed
        update_step_numba(L, spins, neighbors_up, neighbors_down, J, T, h_values, step)
        if step % save_every == 0:
            data.append((step, spins.copy()))
    return data

class Simulation:
    def __init__(self, L=10, T=2.269185, J=1.0, h=None, spins=None):
        self.L = L
        self.T = T
        self.J = J
        self.h = h if h is not None else lambda t, i, j: 0
        self.spins = spins if spins is not None else np.ones((L, L), dtype=int)
        self.neighbors_up, self.neighbors_down = precompute_neighbors_numba(L)

        # Precompute external field values for all steps
        max_steps = 10000  # Adjust this as necessary
        self.h_values = np.zeros((max_steps, L, L))  # Store precomputed fields safely
        for t in range(max_steps):
            for i in range(L):
                for j in range(L):
                    self.h_values[t, i, j] = self.h(t, i, j)

    def simulate(self, steps=100, burn_in_steps=5000, save_every=10):
        """
        Perform the simulation and collect lattice states.
        """
        print(f"[SIMULATE] T={self.T} L={self.L} burn_in={burn_in_steps} steps={steps} save_every={save_every}")
        data = run_simulation(
            self.L,
            steps,
            burn_in_steps,
            save_every,
            self.spins,
            self.neighbors_up,
            self.neighbors_down,
            self.J,
            self.T,
            self.h_values,
        )
        # Prepare DataFrame
        records = [
            {
                "step": step,
                "temperature": self.T,
                "lattice_state_flattened": spins.copy().flatten().tolist(),
                "shape": self.spins.shape,
                "L": self.L,
                "J": self.J,
            }
            for step, spins in data
        ]
        return pd.DataFrame(records)

def main(
    L: int = 100,
    J: float = 1.0,
    steps: int = 100,
    save_every: int = 1,
    output_dir: str = "reports/data/beta",
):
    temperatures = np.linspace(1, 5, 11)
    burn_in_steps = steps // 2  # Half of the total steps
    np.random.seed(42)
    print("[MAIN] Starting simulation")

    df = pd.DataFrame()
    spins = np.random.choice([-1, 1], size=(L, L))  # Random initial configuration

    for temp in temperatures:
        sim = Simulation(L=L, T=temp, J=J, spins=spins.copy())
        result = sim.simulate(steps=steps, burn_in_steps=burn_in_steps, save_every=save_every)
        df = pd.concat([df, result], ignore_index=True)

    os.makedirs(output_dir, exist_ok=True)
    file_path = os.path.join(output_dir, "lattice.json")
    df.to_json(file_path, index=False, orient="records", indent=2)
    print(f"[SAVE] All lattice states saved to {file_path}")

    # Save configuration parameters
    config = {
        "L": L,
        "J": J,
        "steps": steps,
        "burn_in_steps": burn_in_steps,
        "save_every": save_every,
        "temperatures": temperatures.tolist(),
        "output_dir": output_dir
    }
    
    print(config)
    config_path = os.path.join(output_dir, "config.json")
    with open(config_path, "w") as config_file:
        json.dump(config, config_file, indent=2)
    print(f"[SAVE] Configuration saved to {config_path}")

if __name__ == "__main__":
    typer.run(main)
