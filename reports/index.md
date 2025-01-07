# Ising Simulation

```ts
import { reshapeArray, toTidyData } from './components/utils.js'
```

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

```ts
const alpha = await FileAttachment('data/alpha.zip').zip()
```
<!-- 
```ts
display(alpha)
``` -->

```ts
const data = await alpha.file('reports/data/alpha/lattice_states.json').json()
const metrics = await alpha.file('reports/data/alpha/metrics.json').json()
const dynamic_correlations = await alpha.file('reports/data/alpha/dynamic_correlations.json').json()
const spatial_correlations = await alpha.file('reports/data/alpha/spatial_correlations.json').json()
```

<!-- 
```ts
display(data)
``` -->


```ts
const reshapedData = data.map(datum => {
  return {
    ...datum,
    lattice: reshapeArray(datum.lattice_state_flattened, datum.shape[0], datum.shape[1]),
  }
})

const temps = [...new Set(reshapedData.map(d => d.temperature))].reverse()
const aRun = reshapedData.filter(d => d.temperature === temps[0])
```

---

## Lattice state plot

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

## Identifying the critical point

```ts
const metricToPlot = view(Inputs.select(['magnetization', 'average_energy', 'susceptibility', 'specific_heat']))
```

```ts
display(
  Plot.plot({
    marks: [
      Plot.line(metrics, {
        x: "temperature",
        y: metricToPlot
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
  return d.time.map((time, index) => {
      return {
        temperature: d.temperature,
        time: time,
        correlation: d.correlation[index]
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

Dynamic correlations are by comparing each lattice state with the very first state of the simulation (for each temperature). Concretely, the code takes an average of the product ⟨si(0)×si(t)⟩⟨si​(0)×si​(t)⟩ for all lattice spins. This measure gives us a sense of how the initial spin configuration correlates with subsequent configurations over time. If the system is highly ordered and remains that way, the dynamic correlation stays closer to 1. If it quickly randomizes, then the correlation decays toward 0.

When plotting these values, you want to see a correlation curve that starts near 1 (at t=0t=0) and then decays with time. Below the critical temperature, that decay can be slower and decays to a constant value as the system is ordered. Above the critical temperature, we see a faster drop to nearly zero, showing that the system loses memory of its initial arrangement more quickly.

At the critical temperature we see the correlations decay slowly, towards a constant value indicating the existence of long range correlations in the system as is in a state between order and disorder.

### 2. Spatial Correlations

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

Spatial correlations are calculated by shifting the lattice in all possible directions (dx, dy) and measuring how much each site’s spin correlates with its shifted partner. That sum gets binned by the distance ${tex`r=\sqrt{dx^2+dy^2}`}​. In the script, it’s done with a loop that applies np.roll on the lattice to shift it, computes the pairwise product lattice * shifted_lattice, and accumulates results at the correct radial distance. Finally, each bin is normalized by the number of sites that contributed. This gives the correlation value as a function of distance.

When we plot these spatial correlations, you’ll typically see them decay as distance increases. At high temperature, the correlation should fall off quickly, indicating that spins behave more independently. Near or below the critical temperature, the correlation length grows, and we see a slower or even power-law–like decay, indicating large-scale domains of aligned spins.