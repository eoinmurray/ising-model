
# Ising Simulation

```ts
import { reshapeArray, arrayToTidyData } from './components/utils.js';
import { Simulation } from './components/simulation.js';
```

```ts
const data = {
  temperature_energy: await FileAttachment('./data/temperature_energy.json').json(),
  temperature_magnetization: await FileAttachment('./data/temperature_magnetization.json').json(),
  specific_heat: await FileAttachment('./data/specific_heat.json').json(),
  susceptibility: await FileAttachment('./data/susceptibility.json').json(),
  correlations: await FileAttachment('./data/correlations.json').json(),
}
```


<div class="grid grid-cols-2"><div>

```ts
const [temperature, interactionStrength] = [
  view(Inputs.range([0, 5], { value: 1, label: "Temperature", step: 0.1 })),
  view(Inputs.range([0, 5], { value: 1, label: "Interaction Strength", step: 0.1 }))
];
```

```ts
const simulation = new Simulation({ gridSize: 10, temperature, interactionStrength });

const runningSimulation = (function* () {
  for (let i = 0; i<1000; ++i) {
    simulation.update();
    yield {i, simulation};
  }
})();
```

```ts
display(runningSimulation.simulation.getPlots()[0]);
```

</div><div>


- **Lattice**: The system is composed of a lattice of sites (e.g., a 1D line, 2D grid, or 3D cube), where each site holds a spin \( s_i \).
- **Spin States**: Each spin \( s_i \) can take one of two values: \(+1\) (spin up) or \(-1\) (spin down).
- **Interactions**: Neighboring spins interact with each other, with the interaction energy described by the Hamiltonian:

```tex
  H = -J \sum_{\langle i,j \rangle} s_i s_j - h \sum_i s_i
```

- ${tex`J`}: Interaction strength between spins.
- ${tex`\langle i, j \rangle`}: Sum over pairs of neighboring spins.
- ${tex`h`}: External magnetic field applied to the system.

The goal is to study how these spins arrange themselves to minimize energy (at low temperatures) or how they behave statistically at higher temperatures, __in our simulation there is no external field__.

**Description:** Blue (+1) and red (-1) cells show spin configurations. Low temperatures reveal stable clusters, while high temperatures produce random patterns.


</div></div>

---

<div class="grid grid-cols-2"><div>

```ts
display(runningSimulation.simulation.getPlots()[1]);
```

</div><div>

## Hamiltonian Plot

**Key Observations:** Energy decreases sharply initially, stabilizing at equilibrium. High fluctuations suggest transitions.

</div></div>


---

<div class="grid grid-cols-2"><div>


```ts
display(runningSimulation.simulation.getPlots()[2]);
```


</div><div>

## Magnetization Plot

**Key Observations:** At low temperatures, magnetization stabilizes around high values (aligned spins). At high temperatures, it fluctuates near zero due to randomness.

</div></div>

---

<div class="grid grid-cols-2"><div>

```ts
display(Plot.plot({ width: 400, height: 200, marks: [Plot.line(data.temperature_energy, { x: 'temperature', y: 'average_energy' })] }));
```


</div><div>

## Temperature vs. Energy
- **Low Temp:** Spins align, minimizing energy.  
- **High Temp:** Energy plateaus as spins randomize.  
- **Critical Temp:** Sharp rise near Tc indicates phase transition.


</div></div>

---

<div class="grid grid-cols-2"><div>

```ts
display(Plot.plot({ width: 400, height: 200, marks: [Plot.line(data.temperature_magnetization, { x: 'temperature', y: 'average_magnetization' })] }));
```

</div><div>

## Temperature vs. Magnetization

- **Low Temp:** Strong spin alignment yields high |M|.  
- **High Temp:** Randomized spins result in M ≈ 0.  
- **Critical Temp:** Magnetization declines steeply near Tc.

</div></div>

---

<div class="grid grid-cols-2"><div>

```ts
display(Plot.plot({ width: 400, height: 200, marks: [Plot.line(data.specific_heat, { x: 'temperature', y: 'specific_heat' })] }));
```


</div><div>

## Specific Heat
- **Low/High Temp:** Minimal specific heat.  
- **Critical Temp:** Sharp peak due to energy fluctuations.

Specific heat measures how the system's internal energy changes with temperature. It's a thermodynamic property that gives insight into the energy fluctuations of the system. Mathematically, specific heat is computed from the variance of the system's energy.

</div></div>


---

<div class="grid grid-cols-2"><div>

```ts
display(Plot.plot({ width: 400, height: 200, marks: [Plot.line(data.susceptibility, { x: 'temperature', y: 'susceptibility' })] }));
```

</div><div>

## Susceptibility
- **Low/High Temp:** Low susceptibility (ordered/disordered states).  
- **Critical Temp:** Large peak shows enhanced response to external fields.

This measures how sensitive the magnetization is to temperature.

</div></div>

---

<div class="grid grid-cols-2"><div>


```ts
const tidyCorrelations = data.correlations.flatMap(d => 
  d.correlations.map((value, index) => ({
    temperature: d.temperature,
    distance: index + 1,
    correlation: value
  }))
)

display(
  Plot.plot({
    marks: [
      Plot.line(tidyCorrelations, { 
        x: "distance", 
        y: "correlation", 
        stroke: "temperature", 
        strokeWidth: 2, 
        title: "temperature" 
      }),
      Plot.dot(tidyCorrelations, { 
        x: "distance", 
        y: "correlation", 
        fill: "temperature",
        // dx: 5,
        // dy: -5
      })
    ],
    x: { label: "Distance (r)" },
    y: { label: "Correlation C(r)" },
    color: { scheme: "reds", label: "Temperature (T)", legend:true },

    width: 400,
    height: 200
  })
)
```

![](./data/radial_correlations.png)

</div><div>

## Correlations

__Low Temperatures (Below Critical Temperature Tc):__

- High Correlations: Spins are highly ordered (ferromagnetic phase), so the correlation values will start high (close to 1) and decay slowly as the distance increases.
- Long-Range Order: The decay is minimal because spins align over long distances.

__Near the Critical Temperature (T≈Tc​):__

- Slower Decay: Correlations decay more slowly with distance due to critical fluctuations and emerging scale invariance.
- Diverging Correlation Length: Near TcTc​, the system is transitioning from order to disorder, and the correlation length (range over which spins are correlated) grows significantly, possibly spanning the entire lattice.

__High Temperatures (Above Tc):__

- Rapid Decay: Spins become uncorrelated (paramagnetic phase). Correlation values will drop off rapidly to zero as the distance increases.
- Short-Range Order: Correlations exist only over very short distances due to thermal noise dominating spin interactions.

__Expected Features:__

- Curve Shape: At low TT, the curve starts near 1 and decays slowly. At high TT, the curve starts lower and decays rapidly.
- Temperature Dependency: As TT increases, the curves at different temperatures will show progressively faster decay.

</div></div>