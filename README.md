# Ising Model Simulation

This repository contains a Python implementation of the 2D Ising model using the Glauber dynamics. The simulation measures various thermodynamic and statistical properties, such as energy, magnetization, susceptibility, specific heat, and correlations, across different temperatures.

Features:

- Simulate the 2D Ising model on a square lattice.
- Measure thermodynamic properties (e.g., energy, magnetization).
- Compute spatial and temporal correlations.
- Export simulation data in JSON format.
- Plot and save results for visualization using observale plot.



## Tools used

_1. Python_

We run the simulation in python, but instead of using python directly, we use a fancy new tool called [uv](https://docs.astral.sh/uv/) which manages all the dependancies and runs the code in a virtual environment for us. Installation guide below.

_2. Typescript_

We have a live simulation built in typescript running in an Observable Framework notebook. Installation guide below.

## Python Installation

_1. Clone the Repository._

```sh
git clone https://github.com/eoinmurray/ising-model.git
cd ising-model
```

_2. Set Up a Python Environment using [uv](https://docs.astral.sh/uv/getting-started/installation/#installation-methods)._

```sh
brew install uv
```

```sh
uv sync
```

_3. Run the simulation._

```sh
sh run.sh
```

You should see an output like:

```sh
[MAIN] Running simulation
[INIT] Creating Simulation object.
[INIT] Initial spins created.
[SIMULATE] Starting simulation across given temperatures.
  => Setting temperature = 0.1
  => Reinitializing spins for new temperature range.
  => Performing burn-in steps...
  => Collecting data for 1000 steps at T=0.1
  => Setting temperature = 0.35789473684210527
  => Performing burn-in steps...

  ...continued output...
```

- It typically takes the simulation 60 seconds to run.
- It saves all charts and data to the folder `reports/data`.

### Typesrcript and Observable Framework installation

_1. Install nodejs and npm._

Go to the [nodejs website](https://nodejs.org/en) and download the latest version.

_2. Install dependancies._

```sh
npm install
```

_3. Run observable framework_

```
npm run dev
```

The simulation should appear at [http://127.0.0.1:3000/](http://127.0.0.1:3000/)
