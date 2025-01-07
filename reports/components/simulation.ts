import * as Plot from "@observablehq/plot";

type Spin = 1 | -1;
interface LatticePoint {
  x: number;
  y: number;
  spin: Spin;
}

export class Simulation {
  gridSize: number;
  temperature: number;
  interactionStrength: number;
  latticeData: Spin[][];
  isingChart: any;
  hamiltonianHistory: number[];
  magnetizationHistory: number[];
  hamiltonianChart: any;
  magnetizationChart: any;

  constructor({
    gridSize = 10,
    temperature = 1,
    interactionStrength = -1,
  }: {
    gridSize?: number;
    temperature?: number;
    interactionStrength?: number;
  } = {}) {
    this.gridSize = gridSize;
    this.temperature = temperature;
    this.interactionStrength = interactionStrength;
    this.latticeData = this.generateLattice();
    this.hamiltonianHistory = [];
    this.magnetizationHistory = [];

    // Record initial Hamiltonian and magnetization
    this.hamiltonianHistory.push(this.computeHamiltonian());
    this.magnetizationHistory.push(this.computeMagnetization());

    this.isingChart = null;
    this.hamiltonianChart = null;
    this.magnetizationChart = null;
  }

  generateLattice(): Spin[][] {
    return Array.from({ length: this.gridSize }, () =>
      Array.from({ length: this.gridSize }, () => (Math.random() > 0.5 ? 1 : -1))
    );
  }

  computeEnergy(x: number, y: number): number {
    const spin = this.latticeData[y][x];
    const neighbors = [
      this.latticeData[(y - 1 + this.gridSize) % this.gridSize][x],
      this.latticeData[(y + 1) % this.gridSize][x],
      this.latticeData[y][(x - 1 + this.gridSize) % this.gridSize],
      this.latticeData[y][(x + 1) % this.gridSize],
    ];

    if (neighbors.length !== 4) {
      throw new Error("Neighbor calculation error: Unexpected neighbor count");
    }

    const neighborSum = neighbors.reduce((sum, neighborSpin) => sum + neighborSpin, 0);
    // Local energy for this spin site = -J * s_i * sum_{neighbors}(s_j)
    return -this.interactionStrength * spin * neighborSum;
  }

  computeHamiltonian(): number {
    let totalEnergy = 0;
    for (let y = 0; y < this.gridSize; y++) {
      for (let x = 0; x < this.gridSize; x++) {
        totalEnergy += this.computeEnergy(x, y);
      }
    }
    return totalEnergy / 2; // Avoid double counting each bond
  }

  computeMagnetization(): number {
    return this.latticeData.flat().reduce((sum, spin) => sum + spin, 0);
  }

  /**
   * One Monte Carlo step (Glauber dynamics).
   * We pick a random spin and flip it with probability:
   *    p_flip = 1 / (1 + exp( ΔE / T ))
   */
  update(): void {
    // Pick a random site
    const x = Math.floor(Math.random() * this.gridSize);
    const y = Math.floor(Math.random() * this.gridSize);

    const spin = this.latticeData[y][x];
    // Sum of neighbor spins
    const neighbors = [
      this.latticeData[(y - 1 + this.gridSize) % this.gridSize][x],
      this.latticeData[(y + 1) % this.gridSize][x],
      this.latticeData[y][(x - 1 + this.gridSize) % this.gridSize],
      this.latticeData[y][(x + 1) % this.gridSize],
    ];
    const neighborSum = neighbors.reduce((sum, neighborSpin) => sum + neighborSpin, 0);

    // Current local energy = -J * spin * neighborSum
    const energyBeforeFlip = -this.interactionStrength * spin * neighborSum;
    // After flipping spin, local energy = -J * (-spin) * neighborSum = +J * spin * neighborSum
    const energyAfterFlip = -energyBeforeFlip;
    // ∆E = E_after - E_before = -energyBeforeFlip - energyBeforeFlip = -2 * energyBeforeFlip
    const energyChange = energyAfterFlip - energyBeforeFlip;

    // Glauber acceptance probability
    const flipProbability = 1 / (1 + Math.exp(energyChange / this.temperature));
    if (Math.random() < flipProbability) {
      this.latticeData[y][x] *= -1;
    }

    // Record Hamiltonian and magnetization
    this.hamiltonianHistory.push(this.computeHamiltonian());
    this.magnetizationHistory.push(this.computeMagnetization());
  }

  createPlot(): any {
    const data: LatticePoint[] = [];
    for (let y = 0; y < this.gridSize; y++) {
      for (let x = 0; x < this.gridSize; x++) {
        data.push({ x, y, spin: this.latticeData[y][x] });
      }
    }

    return Plot.plot({
      marks: [
        Plot.cell(data, {
          x: "x",
          y: "y",
          fill: (d: LatticePoint) => (d.spin === 1 ? "blue" : "red"),
        }),
      ],
      width: 400,
      height: 400,
    });
  }

  createHamiltonianPlot(): any {
    return Plot.plot({
      marks: [
        Plot.line(
          this.hamiltonianHistory.map((value, index) => ({ index, value })),
          { x: "index", y: "value" }
        ),
      ],
      width: 400,
      height: 200,
      x: { label: "Time Step" },
      y: { label: "Hamiltonian" },
    });
  }

  createMagnetizationPlot(): any {
    return Plot.plot({
      marks: [
        Plot.line(
          this.magnetizationHistory.map((value, index) => ({ index, value })),
          { x: "index", y: "value" }
        ),
      ],
      width: 400,
      height: 200,
      x: { label: "Time Step" },
      y: { label: "Magnetization" },
    });
  }

  /**
   * Generate the three plots: the spins (blue=+1, red=-1),
   * Hamiltonian over time, and magnetization over time.
   */
  getPlots(): any[] {
    this.isingChart = this.createPlot();
    this.hamiltonianChart = this.createHamiltonianPlot();
    this.magnetizationChart = this.createMagnetizationPlot();

    return [this.isingChart, this.hamiltonianChart, this.magnetizationChart];
  }
}
