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

```ts
display(alpha)
```

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

## Identifying the critical point

We can plot the following graphs: magnetization, average_energy, susceptibility, and specific_heat.
- For magnetization we should see magnetization drop at the critical point
- For average_energy we should see the center of an ${tex`S`} curve at the critical point.
- For susceptibility and specific_heat we should see a peak at the critical point.

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