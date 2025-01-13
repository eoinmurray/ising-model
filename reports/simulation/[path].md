# Ising Simulation Run: ${observable.params.path}

> _Note the plots may take 30s to load as the data is loaded._

```ts
import { reshapeArray, toTidyData } from '../components/utils.js'
```

```ts
const div = html`<div style="background:red;width:0;height:13px;">`;
(async () => {
  let sum = 0;
  for (let i = 0, n = 500; i < n; ++i) {
    sum += i * i;
    div.style.width = (i / n) * 100 + "%";
    await new Promise(requestAnimationFrame);
  }
  div.value = sum;
  // div.dispatchEvent(new CustomEvent("input"));
  // div.style.background = "black";
})();

display(div)
```

```ts
display("Data is loading.")
const data = await FileAttachment(`../data/${observable.params.path}/data.zip`).zip()
display(data)
const config = await data.file('config.json').json()
const lattice = await data.file('lattice.json').json()
const metrics = await data.file('metrics.json').json()
const dynamic_correlations = await data.file('dynamic_autocorrelations.json').json()
const correlation_times = await data.file('correlation_times.json').json()
const spatial_correlations = await data.file('spatial_correlations.json').json()
display("Data loaded.")
div.style.width = 100 + "%";
div.style.display = "none";
display(config)
```

The Ising is a tool used in statistical physics to study phase transitions and collective dynamics in systems with
many interacting parts. Initially designed to represent ferromagnetism, it models spins on a lattice that can take one of two states, interacting with their neighbors. Its applications extend to areas like neural network activity, protein folding and quantum cluster states. The model's simplicity allows for analytical solutions in certain cases and computational approaches in others, making it broadly applicable for exploring the connection between microscopic interactions and macroscopic phenomena.

Here I'm making a small monte carlo simulation of an Ising Model with Glauber dynamics. 

### 1. **Energy of the Ising Model**:
The energy of a configuration of spins ${tex`\{ \sigma_i \}`} (where ${tex`\sigma_i = \pm 1`}) on a lattice is given by:

```tex
E = -J \sum_{\langle i,j \rangle} \sigma_i \sigma_j - h \sum_i \sigma_i
```

Where:
- ${tex` J `} is the interaction strength between neighboring spins (positive for ferromagnetic, negative for antiferromagnetic),
- ${tex` \langle i,j \rangle `} denotes summation over nearest neighbors,
- ${tex` h `} is an external magnetic field, and
- ${tex` \sigma_i `} is the spin at site ${tex` i `}.

### 2. **Glauber Dynamics Update Rule**:
The Glauber dynamics updates the spin at site ${tex` i `} based on the local energy configuration. The probability of a spin ${tex` \sigma_i `} flipping to ${tex` \sigma_i' `} is:

```tex
P(\sigma_i \to \sigma_i') = \frac{1}{1 + \exp\left( \frac{\Delta E}{k_B T} \right)}
```

Where:
- ${tex` \Delta E = E(\sigma_i' \text{ flipped}) - E(\sigma_i) `} is the change in energy when the spin flips,
- ${tex` k_B `} is the Boltzmann constant,
- ${tex` T `} is the temperature.

For the case of a single spin flip ${tex` \sigma_i \to \sigma_i' = -\sigma_i `}, the energy change ${tex` \Delta E `} for the flip is:

```tex
\Delta E = 2 \sigma_i \left( \sum_{\langle i,j \rangle} J \sigma_j + h \right)
```

This means that the probability of flipping a spin depends on the interaction with neighboring spins and the external field. At high temperatures (${tex` T \to \infty `}), flips occur more freely, while at low temperatures (${tex` T \to 0 `}), the system tends to settle into aligned configurations.

Note the critical temperature for an ising model should lie around 

```tex
T_c = \frac{2J}{ln(1 + \sqrt{2})} \approx 2.27 J
```

in our sim J is always 1.0, so the critical temp is ${tex`\approx`} 2.27.

---

## Lattice state plot

```ts
const reshapedData = lattice.map(datum => {return {...datum,lattice: reshapeArray(datum.lattice_state_flattened, datum.shape[0], datum.shape[1]),}})
const temps = [...new Set(reshapedData.map(d => d.temperature))].reverse()
const aRun = reshapedData.filter(d => d.temperature === temps[0])
```

```ts
const selectedTemp = view(Inputs.select(temps, { label: "Temperature" }));
```

```ts
const run = reshapedData.filter(d => d.temperature === selectedTemp)
const stepSize = run[1].step - run[0].step
```

```ts
const step = view(Inputs.range([0, (aRun.length-1)*stepSize], { step: stepSize, value: 0, label: "Time Step Index" }))
```

```ts
display(
  Plot.plot({
    width: 400,
    height: 400,
    padding: 0,
    marks: [
      Plot.cell(toTidyData(run[step / stepSize].lattice), {
        x: 'row',
        y: 'column',
        fill: 'value'
      })
    ],
    color: {
      domain: [0, 1]
    },
    x: { axis: null },
    y: { axis: null },
  })
)
```

Here we plot the lattice state for various temperatures and time steps, play around with the sliders to 
view the dynamics behaviour of the lattice.

The simulation calculates the lattice states using Glauber dynamics in the `update_step_numba` function. At each step, a random site is picked, and the change in energy (${tex`\delta E`}) due to flipping that spin is computed from the local interactions (the sum of the nearest neighbors’ spins and any external field). The flip is then accepted with probability ${tex`p_{flip} = \frac{1}{1+e^{\delta E / T}}`}​. Repeating this process for many steps yields a time evolution of the spin configuration.

At low temperatures, we often see large, stable clusters where spins align and rarely flip, indicating an ordered phase. Near critical temperature (${tex`T_c`}), these clusters form and dissolve in a more balanced manner, and we catch glimpses of interesting domain structures. At high temperatures, spins will flip frequently and produce a more disordered or “noisy” lattice with low magnetization.

---

## Correlation times

```ts
display(
  Plot.plot({
    marks: [
      Plot.line(correlation_times, {
        x: "temperature",
        y: "char_time"
      })
    ]
  })
)

```

## Identifying the critical point


```ts
display(
  Plot.plot({
    marks: [
      Plot.line(metrics, {
        x: "temperature",
        y: "magnetization"
      })
    ]
  })
)
```


```ts
display(
  Plot.plot({
    marks: [
      Plot.line(metrics, {
        x: "temperature",
        y: "average_energy"
      })
    ]
  })
)
```


```ts
display(
  Plot.plot({
    marks: [
      Plot.line(metrics, {
        x: "temperature",
        y: "susceptibility"
      })
    ]
  })
)
```


```ts
display(
  Plot.plot({
    marks: [
      Plot.line(metrics, {
        x: "temperature",
        y: "specific_heat"
      })
    ]
  })
)
```

We can plot the following graphs: magnetization, average_energy, susceptibility, and specific_heat.
- For magnetization we should see magnetization drop at the critical point
- For average_energy we should see the center of an ${tex`S`} curve at the critical point.
- For susceptibility and specific_heat we should see a peak at the critical point.

---

## Correlations

```ts
const tidyDynamicCorrelations = dynamic_correlations.flatMap(d => {
  return d.time_lag.map((time, index) => {
      return {
        temperature: d.temperature,
        time,
        correlation: d.autocorrelation[index]
      }
    })
})

const tidySpatialCorrelations = spatial_correlations.flatMap(d => {
  return d.distance.map((distance, index) => {
      return {
        temperature: d.temperature,
        distance: distance,
        correlation: d.correlation[index]
      }
    })
})
```


### 1. Dynamic Correlations



```ts
display(
  Plot.plot({
    marks: [
      Plot.lineY(tidyDynamicCorrelations, {
        x: 'time',
        y: 'correlation',
        stroke: 'temperature',
        tip: 'x'
      })
    ],
    x: {
      domain: [0, 25]
    },
    y: {
      domain: [0, 1]
    },
    color: {
      scheme: 'reds',
      legend: true
    }
  })
)
```

Dynamic correlations are by comparing each lattice state with the very first state of the simulation (for each temperature). The code takes an average of the product ${tex`\langle s_i(0) \times s_i(t) \rangle`} for all lattice spins. This measure gives us a sense of how the initial spin configuration correlates with subsequent configurations over time. If the system is highly ordered and remains that way, the dynamic correlation stays closer to 1. If it quickly randomizes, then the correlation decays toward 0.

When plotting these values, you want to see a correlation curve that starts near 1 (at t=0t=0) and then decays with time. Below the critical temperature, that decay can be slower and decays to a constant value as the system is ordered. Above the critical temperature, we see a faster drop to nearly zero, showing that the system loses memory of its initial arrangement more quickly.

At the critical temperature we see the correlations decay slowly, towards a constant value indicating the existence of long range correlations in the system as is in a state between order and disorder.

```ts
display(
  Plot.plot({
    marks: [
      Plot.lineY(tidySpatialCorrelations, {
        x: 'distance',
        y: 'correlation',
        stroke: 'temperature',
        tip: 'x'
      })
    ],
    x: {
      domain: [0, 50]
    },
    y: {
      domain: [0, 1]
    },
    color: {
      scheme: 'reds',
      legend: true
    }
  })
)
```

Spatial correlations are calculated by shifting the last step lattice in all possible directions (dx, dy) and measuring how much each site’s spin correlates with its shifted partner. That sum gets binned by the distance ${tex`r=\sqrt{dx^2+dy^2}`}​. In the script, it’s done with a loop that applies np.roll on the lattice to shift it, computes the pairwise product lattice * shifted_lattice, and accumulates results at the correct radial distance. Finally, each bin is normalized by the number of sites that contributed. This gives the correlation value as a function of distance.

When we plot these spatial correlations, you’ll typically see them decay as distance increases. At high temperature, the correlation should fall off quickly, indicating that spins behave more independently. Near or below the critical temperature, the correlation length grows, and we see a slower or even power-law–like decay, indicating large-scale domains of aligned spins.