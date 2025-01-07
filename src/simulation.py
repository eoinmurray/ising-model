import faulthandler
faulthandler.enable()

import numpy as np
import pandas as pd
import time
import os
from numba import njit

@njit
def precompute_neighbors_numba(L):
    indices = np.arange(L)
    neighbors_up = np.roll(indices, -1)
    neighbors_down = np.roll(indices, 1)
    return neighbors_up, neighbors_down

@njit
def update_step_numba(L, spins, neighbors_up, neighbors_down, J, T, h_values, step):
    """
    Perform a single update step using Glauber dynamics.
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
        dE = 2 * J * spins[i, j] * (neighbors_sum + field)
        p_flip = 1.0 / (1.0 + np.exp(dE / T))

        if np.random.rand() < p_flip:
            spins[i, j] *= -1

@njit
def run_simulation(L, steps, save_every, spins, neighbors_up, neighbors_down, J, T, h_values):
    """
    Run the simulation for the specified number of steps.
    """
    data = []
    for step in range(steps):
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
        max_steps = 1000  # Adjust this as necessary
        self.h_values = np.zeros((max_steps, L, L))  # Store precomputed fields safely
        for t in range(max_steps):
            for i in range(L):
                for j in range(L):
                    self.h_values[t, i, j] = self.h(t, i, j)

    def simulate(self, steps=100, save_every=10):
        """
        Perform the simulation and collect lattice states.
        """
        print(f"[SIMULATE] Running simulation for T={self.T} steps={steps} save_every={save_every}")
        data = run_simulation(
            self.L,
            steps,
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
                "lattice_state_flattened": spins.flatten().tolist(),
                "shape": self.spins.shape,
                "L": self.L,
                "J": self.J,
            }
            for step, spins in data
        ]
        return pd.DataFrame(records)

def encoded_signal(time, i, j):
    """
    Example of an external oscillating field.
    """
    center_x, center_y = 50, 50  # Center of the signal
    radius = 10                # Affected region
    strength = 10.0            # Field strength

    distance = np.sqrt((i - center_x)**2 + (j - center_y)**2)
    if distance <= radius:
        return strength * np.sin(2 * np.pi * time / 50.0)  # Oscillating field
    return 0

def metrics(input_file: str = "reports/data/alpha/lattice_states.json"):
    """
    Compute metrics from the simulation output.

    Args:
        input_file: Path to the lattice states file.
    """
    print("[METRICS] Computing metrics")
    start_time = time.time()

    df = pd.read_json(input_file)
    metrics = []

    max_step = df['step'].unique().max()
    df_max_step = df[df['step'] == max_step]
    for _, row in df_max_step.iterrows():
        lattice = np.array(row["lattice_state_flattened"]).reshape(row["shape"])
        J = row["J"]

        magnetization = np.abs(np.sum(lattice)) / lattice.size

        energy = -J * np.sum(lattice * np.roll(lattice, 1, axis=0)) \
                 -J * np.sum(lattice * np.roll(lattice, 1, axis=1))
        average_energy = energy / lattice.size
        specific_heat = (np.var(lattice) * J**2) / (row["temperature"]**2 * lattice.size)

        susceptibility = (np.var(lattice) * J) / (row["temperature"] * lattice.size)

        metrics.append({
            "temperature": row["temperature"],
            "magnetization": magnetization,
            "average_energy": average_energy,
            "specific_heat": specific_heat,
            "susceptibility": susceptibility
        })

    metrics_df = pd.DataFrame(metrics)

    output_file = input_file.replace("lattice_states.json", "metrics.json")
    metrics_df.to_json(output_file, index=False, orient="records", indent=2)
    print(metrics_df)
    print(f"[SAVE] All metrics saved to {output_file}")

    print("[METRICS] Computation complete.")
    end_time = time.time()
    print(f"[MAIN] Execution time: {end_time - start_time:.6f} seconds")

def spatial_correlations(input_file: str = "reports/data/alpha/lattice_states.json"):
  """
  Calculate spatial correlations from the lattice states.

  Args:
    input_file: Path to the lattice states file.
  """
  print("[SPATIAL_CORRELATION] Calculating spatial correlations")
  start_time = time.time()

  def calculate_correlation(lattice, L):
    print("[CALCULATE_CORRELATION] Optimized correlation calculation")
    max_distance = int(np.ceil(np.sqrt(2) * L))
    correlation = np.zeros(max_distance)
    count = np.zeros(max_distance)

    for dx in range(-L + 1, L):
      for dy in range(-L + 1, L):
        r = int(round(np.sqrt(dx**2 + dy**2)))
        if 0 <= r < max_distance:
          shifted_lattice = np.roll(np.roll(lattice, dx, axis=0), dy, axis=1)
          correlation[r] += np.sum(lattice * shifted_lattice)
          count[r] += lattice.size

    correlation /= count
    return correlation

  df = pd.read_json(input_file)
  correlations = []

  for temperature in df['temperature'].unique():
    print(f"[SPATIAL_CORRELATION] Processing temperature {temperature}")
    temp_df = df[df['temperature'] == temperature]
    max_step = temp_df['step'].max()
    last_state = temp_df[temp_df['step'] == max_step].iloc[0]

    lattice = np.array(last_state["lattice_state_flattened"]).reshape(last_state["shape"])
    print("[SPATIAL_CORRELATION] Lattice reshaped for analysis")
    L = last_state["L"]
    corr = calculate_correlation(lattice, L)

    correlations.append({
      "temperature": temperature,
      "distance": np.arange(len(corr)).tolist(),
      "correlation": corr.tolist()
    })

  correlation_df = pd.DataFrame(correlations)

  output_file = input_file.replace("lattice_states.json", "spatial_correlations.json")
  os.makedirs(os.path.dirname(output_file), exist_ok=True)
  correlation_df.to_json(output_file, index=False, orient="records", indent=2)

  print("[SPATIAL_CORRELATION] Correlation data saved to DataFrame and JSON")
  print(correlation_df)
  print(f"[SAVE] Spatial correlations saved to {output_file}")

  print("[SPATIAL_CORRELATION] Calculation complete.")
  end_time = time.time()
  print(f"[MAIN] Execution time: {end_time - start_time:.6f} seconds")

def dynamic_correlations(input_file: str = "reports/data/alpha/lattice_states.json"):
  """
  Calculate dynamic correlations across steps for each temperature.

  Args:
    input_file: Path to the lattice states file.
  """
  print("[DYNAMIC_CORRELATION] Calculating dynamic correlations")
  start_time = time.time()

  def calculate_dynamic_correlation(lattice_states):
    print("[CALCULATE_DYNAMIC] Computing time-based correlations")
    num_steps = len(lattice_states)
    correlations = np.zeros(num_steps)

    # We compare all states with the first, as a reference
    for t in range(num_steps):
      correlations[t] = np.mean(lattice_states[0] * lattice_states[t])

    return correlations

  df = pd.read_json(input_file)
  dynamic_correlations = []

  for temperature in df['temperature'].unique():
    print(f"[DYNAMIC_CORRELATION] Processing temperature {temperature}")
    temp_df = df[df['temperature'] == temperature]
    lattice_states = [
      np.array(row["lattice_state_flattened"]).reshape(row["shape"])
      for _, row in temp_df.iterrows()
    ]

    correlations = calculate_dynamic_correlation(lattice_states)

    dynamic_correlations.append({
      "temperature": temperature,
      "time": np.arange(len(correlations)).tolist(),
      "correlation": correlations.tolist()
    })

  dynamic_correlation_df = pd.DataFrame(dynamic_correlations)

  output_file = input_file.replace("lattice_states.json", "dynamic_correlations.json")
  os.makedirs(os.path.dirname(output_file), exist_ok=True)
  dynamic_correlation_df.to_json(output_file, index=False, orient="records", indent=2)

  print("[DYNAMIC_CORRELATION] Correlation data saved to DataFrame and JSON")
  print(dynamic_correlation_df)
  print(f"[SAVE] Dynamic correlations saved to {output_file}")

  print("[DYNAMIC_CORRELATION] Calculation complete.")
  end_time = time.time()
  print(f"[MAIN] Execution time: {end_time - start_time:.6f} seconds")


def main():
    np.random.seed(42)
    temperatures = [1, 1.5, 2, 2.27, 2.5, 3, 3.5, 4, 4.5, 5]
    L=100
    J=1.0
    steps=100000
    save_every= int(steps / 1000)
    output_dir="reports/data/alpha"

    print("[MAIN] Starting simulation")

    df = pd.DataFrame()
    
    # spins = -np.ones((L, L))  # Start with all spins down
    spins = np.random.choice([-1, 1], size=(L, L))  # Initialize spins randomly

    for temp in temperatures:
        result = Simulation(
            L=L, 
            T=temp, 
            J=J, 
            spins=spins.copy(),
            # h=encoded_signal
        ).simulate(steps=steps, save_every=save_every)
        df = pd.concat([df, result], ignore_index=True)

    print(df)
    os.makedirs(output_dir, exist_ok=True)
    file_path = os.path.join(output_dir, "lattice_states.json")
    
    df.to_json(file_path, index=False, orient="records", indent=2)

    print(df)
    print(f"[SAVE] All lattice states saved to {file_path}")

    # Compute metrics
    metrics(input_file="reports/data/alpha/lattice_states.json")

    # Calculate spatial correlations
    spatial_correlations(input_file="reports/data/alpha/lattice_states.json")

    # Calculate dynamic correlations
    dynamic_correlations(input_file="reports/data/alpha/lattice_states.json")
    
if __name__ == "__main__":
    main()
