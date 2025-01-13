import faulthandler
faulthandler.enable()

import numpy as np
import pandas as pd
import time
import os
import typer
from numba import njit, prange

# For plotting
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

sns.set(style="whitegrid")


def metrics(input_file: str):
    """
    Compute metrics from the simulation output and produce summary plots.

    Args:
        input_file: Path to the lattice states file.
    """
    print("[METRICS] Computing metrics")
    start_time = time.time()

    df = pd.read_json(input_file)
    metrics_data = []

    max_step = df['step'].unique().max()
    df_max_step = df[df['step'] == max_step]

    for _, row in df_max_step.iterrows():
        lattice = np.array(row["lattice_state_flattened"]).reshape(row["shape"])
        J = row["J"]

        # Calculate magnetization
        magnetization = np.abs(np.sum(lattice)) / lattice.size

        # Calculate average energy
        energy = -J * np.sum(lattice * np.roll(lattice, 1, axis=0)) \
                 -J * np.sum(lattice * np.roll(lattice, 1, axis=1))
        average_energy = energy / lattice.size

        # Calculate specific heat
        specific_heat = (np.var(lattice) * J**2) / (row["temperature"]**2 * lattice.size)

        # Calculate susceptibility
        susceptibility = (np.var(lattice) * J) / (row["temperature"] * lattice.size)

        metrics_data.append({
            "temperature": row["temperature"],
            "magnetization": magnetization,
            "average_energy": average_energy,
            "specific_heat": specific_heat,
            "susceptibility": susceptibility
        })

    metrics_df = pd.DataFrame(metrics_data)

    # Save JSON
    output_file = input_file.replace("lattice.json", "metrics.json")
    metrics_df.to_json(output_file, index=False, orient="records", indent=2)
    print(metrics_df)
    print(f"[SAVE] All metrics saved to {output_file}")

    # -----------------
    # Create plots
    # -----------------
    # 1. Magnetization vs Temperature
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.lineplot(
        data=metrics_df, x="temperature", y="magnetization", marker="o", ax=ax
    )
    ax.set_title("Magnetization vs Temperature")
    ax.set_xlabel("Temperature")
    ax.set_ylabel("Magnetization")
    fig.tight_layout()
    plot_file = output_file.replace(".json", "_magnetization.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()

    # 2. Average Energy vs Temperature
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.lineplot(
        data=metrics_df, x="temperature", y="average_energy", marker="o", ax=ax
    )
    ax.set_title("Average Energy vs Temperature")
    ax.set_xlabel("Temperature")
    ax.set_ylabel("Average Energy")
    fig.tight_layout()
    plot_file = output_file.replace(".json", "_avg_energy.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()

    # 3. Specific Heat vs Temperature
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.lineplot(
        data=metrics_df, x="temperature", y="specific_heat", marker="o", ax=ax
    )
    ax.set_title("Specific Heat vs Temperature")
    ax.set_xlabel("Temperature")
    ax.set_ylabel("Specific Heat")
    fig.tight_layout()
    plot_file = output_file.replace(".json", "_specific_heat.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()

    # 4. Susceptibility vs Temperature
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.lineplot(
        data=metrics_df, x="temperature", y="susceptibility", marker="o", ax=ax
    )
    ax.set_title("Susceptibility vs Temperature")
    ax.set_xlabel("Temperature")
    ax.set_ylabel("Susceptibility")
    fig.tight_layout()
    plot_file = output_file.replace(".json", "_susceptibility.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()

    print("[METRICS] Computation and plotting complete.")
    end_time = time.time()
    print(f"[MAIN] Execution time: {end_time - start_time:.6f} seconds")


def spatial_correlations(input_file: str):
    """
    Calculate spatial correlations from the lattice states and plot them.

    Args:
        input_file: Path to the lattice states file.
    """
    print("[SPATIAL_CORRELATION] Calculating spatial correlations")
    start_time = time.time()

    def calculate_correlation(lattice, L):
        """
        A direct calculation of correlation vs distance.
        """
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

    # We might gather all temperatures and store correlation vs distance
    for temperature in df['temperature'].unique():
        print(f"[SPATIAL_CORRELATION] Processing temperature {temperature}")
        temp_df = df[df['temperature'] == temperature]
        max_step = temp_df['step'].max()
        last_state = temp_df[temp_df['step'] == max_step].iloc[0]

        lattice = np.array(last_state["lattice_state_flattened"]).reshape(last_state["shape"])
        L = last_state["L"]
        corr = calculate_correlation(lattice, L)

        correlations.append({
            "temperature": temperature,
            "distance": np.arange(len(corr)).tolist(),
            "correlation": corr.tolist()
        })

    correlation_df = pd.DataFrame(correlations)

    # Save to JSON
    output_file = input_file.replace("lattice.json", "spatial_correlations.json")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    correlation_df.to_json(output_file, index=False, orient="records", indent=2)

    print("[SPATIAL_CORRELATION] Correlation data saved to DataFrame and JSON")
    print(correlation_df)
    print(f"[SAVE] Spatial correlations saved to {output_file}")

    # -----------------
    # Create plots
    # -----------------
    # We'll plot correlation vs distance for each temperature on the same figure
    fig, ax = plt.subplots(figsize=(7, 5))

    # Iterate over the correlation data
    for idx, row in correlation_df.iterrows():
        distances = row["distance"]
        corr_vals = row["correlation"]
        temp = row["temperature"]
        ax.plot(distances, corr_vals, label=f"T={temp}")

    ax.set_title("Spatial Correlation vs Distance")
    ax.set_xlabel("Distance (r)")
    ax.set_ylabel("Correlation")
    ax.legend(loc="best")
    fig.tight_layout()

    plot_file = output_file.replace(".json", "_plot.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()

    print("[SPATIAL_CORRELATION] Calculation and plotting complete.")
    end_time = time.time()
    print(f"[MAIN] Execution time: {end_time - start_time:.6f} seconds")


@njit(parallel=True)
def compute_autocorrelation(lattice_states_flat, L, num_states, max_tau):
    """
    Compute the autocorrelation function up to 'max_tau' for a series of flattened lattices
    across steps.

    Args:
        lattice_states_flat: 2D array of shape [num_states, #spins_per_lattice].
        L: Lattice dimension (not used here, but left for completeness).
        num_states: Total number of states in time dimension.
        max_tau: Maximum time lag for the correlation calculation.

    Returns:
        1D array of autocorrelations (size = max_tau).
    """
    autocorr = np.zeros(max_tau)
    mean_mag = np.mean(lattice_states_flat)

    for tau in prange(max_tau):
        sum_corr = 0.0
        count = 0
        for t in prange(num_states - tau):
            dev_t = lattice_states_flat[t] - mean_mag
            dev_t_tau = lattice_states_flat[t + tau] - mean_mag
            sum_corr += np.mean(dev_t * dev_t_tau)
            count += 1
            
        if count > 0:
            autocorr[tau] = sum_corr / count
            
    # Normalize by tau=0 value
    if autocorr[0] != 0:
        autocorr = autocorr / autocorr[0]
    return autocorr


def _exponential_decay(t, tau_c):
    """
    Single-parameter exponential decay function for fitting:
        C(t) = exp(-t / tau_c).
    """
    return np.exp(-t / tau_c)


def _fit_correlation_time(time_lag, autocorr):
    """
    Estimate the characteristic correlation time tau_c by fitting the
    autocorrelation data to a simple exponential decay, ignoring any offset.
    """
    # Convert to numpy arrays
    time_lag_array = np.array(time_lag, dtype=float)
    autocorr_array = np.array(autocorr, dtype=float)

    # We skip the first point (tau=0 => autocorr=1) to avoid log(1) = 0 dominating the fit
    # Also remove any negative or zero points that might come from numerical issues
    valid_mask = (time_lag_array > 0) & (autocorr_array > 0)
    if not np.any(valid_mask):
        # If no valid points, return a default
        return 0.0

    time_lag_array = time_lag_array[valid_mask]
    autocorr_array = autocorr_array[valid_mask]

    # Perform a curve fit to C(t) = exp(-t / tau_c)
    try:
        popt, _ = curve_fit(_exponential_decay, time_lag_array, autocorr_array, p0=[1.0])
        tau_c = popt[0]
    except RuntimeError:
        # If the fit fails, return a fallback
        tau_c = 0.0

    return tau_c


def dynamic_correlations(input_file: str, max_tau_ratio=0.5):
    """
    Calculate dynamic autocorrelations across steps for each temperature, 
    estimate the characteristic correlation time, and produce plots.

    Args:
        input_file: Path to the lattice states file.
        max_tau_ratio: Ratio of max_tau to total number of states (e.g., 0.1 for 10%).
    """
    print("[DYNAMIC_CORRELATION] Calculating dynamic autocorrelations")
    start_time = time.time()

    # Load the lattice states
    df = pd.read_json(input_file)
    dynamic_correlations_data = []

    # We'll store [temperature, char_time] to plot correlation time vs temperature
    correlation_times = []

    # Iterate over each unique temperature
    for temperature in df['temperature'].unique():
        print(f"[DYNAMIC_CORRELATION] Processing temperature {temperature}")
        temp_df = df[df['temperature'] == temperature].sort_values('step')
        num_states = len(temp_df)
        if num_states < 2:
            print(f"[WARNING] Not enough states for temperature {temperature} to compute correlations.")
            continue

        # Extract and flatten all lattice states
        lattice_states = temp_df["lattice_state_flattened"].tolist()
        lattice_states_flat = np.array([np.array(state).flatten() for state in lattice_states])

        # Determine the maximum tau based on the ratio
        max_tau = int(num_states * max_tau_ratio)
        max_tau = max_tau if max_tau > 0 else 1  # Ensure at least 1

        # Compute autocorrelation using the optimized function
        autocorr = compute_autocorrelation(lattice_states_flat,
                                           temp_df.iloc[0]["L"],
                                           num_states, 
                                           max_tau)

        # Fit or estimate characteristic correlation time (tau_c)
        time_lags = np.arange(max_tau)
        tau_c = _fit_correlation_time(time_lags, autocorr)

        correlation_times.append({"temperature": temperature, "char_time": tau_c})

        dynamic_correlations_data.append({
            "temperature": temperature,
            "time_lag": time_lags.tolist(),
            "autocorrelation": autocorr.tolist(),
            "char_time": tau_c
        })

    # Create a DataFrame from the results
    dynamic_correlation_df = pd.DataFrame(dynamic_correlations_data)

    # Save to JSON
    output_file = input_file.replace("lattice.json", "dynamic_autocorrelations.json")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    dynamic_correlation_df.to_json(output_file, index=False, orient="records", indent=2)

    print("[DYNAMIC_CORRELATION] Autocorrelation data saved to DataFrame and JSON")
    print(dynamic_correlation_df)
    print(f"[SAVE] Dynamic autocorrelations saved to {output_file}")

    # -----------------
    # Create plots
    # -----------------

    # 1) Plot autocorrelation vs time_lag for each temperature
    fig, ax = plt.subplots(figsize=(7, 5))
    for idx, row in dynamic_correlation_df.iterrows():
        time_lag = row["time_lag"]
        corr_vals = row["autocorrelation"]
        temp = row["temperature"]
        ax.plot(time_lag, corr_vals, marker="o", label=f"T={temp} (tau_c={row['char_time']:.2f})")

    ax.set_title("Autocorrelation vs Time Lag")
    ax.set_xlabel("Time Lag (tau)")
    ax.set_ylabel("Autocorrelation")
    ax.legend(loc="best")
    fig.tight_layout()

    plot_file = output_file.replace(".json", "_plot.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()

    # 2) Plot characteristic correlation time vs temperature
    correlation_times_df = pd.DataFrame(correlation_times)
    
    # Save correlation times to JSON
    correlation_times_output_file = input_file.replace("lattice.json", "correlation_times.json")
    correlation_times_df.to_json(correlation_times_output_file, index=False, orient="records", indent=2)
    print(f"[SAVE] Correlation times saved to {correlation_times_output_file}")
    
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.scatterplot(
        data=correlation_times_df, 
        x="temperature", 
        y="char_time", 
        ax=ax, 
        s=60
    )
    sns.lineplot(
        data=correlation_times_df.sort_values("temperature"), 
        x="temperature", 
        y="char_time", 
        ax=ax
    )
    ax.set_title("Characteristic Correlation Time vs Temperature")
    ax.set_xlabel("Temperature")
    ax.set_ylabel("Correlation Time (tau_c)")
    fig.tight_layout()

    plot_file = output_file.replace(".json", "_char_time_plot.png")
    plt.savefig(plot_file, dpi=300)
    plt.close()

    print("[DYNAMIC_CORRELATION] Calculation and plotting complete.")
    end_time = time.time()
    print(f"[MAIN] Execution time: {end_time - start_time:.6f} seconds")


def main(input_file = "reports/data/beta/lattice.json"):
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(input_file), exist_ok=True)

    metrics(input_file)
    spatial_correlations(input_file)
    dynamic_correlations(input_file)


if __name__ == "__main__":
    typer.run(main)
