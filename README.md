# Photonic Ring Resonator ODE Solver and Monte Carlo Analysis

This repository provides MATLAB scripts for simulating high-order ordinary differential equation (ODE) solvers implemented via cascaded microring resonators (MRRs) and performing Monte Carlo analysis to assess performance under parameter variations.

## 🗂️ Repository Structure

```
.
├── ODE_solver.m       % Main script: sets up parameters, computes ODE solutions, and plots results
├── MonteCarlo.m       % Function: runs Monte Carlo trials to analyze perturbations in coupling, loss, and detuning
└── README.md          % Project overview and usage instructions
```

## 📋 Prerequisites

* MATLAB R2018a or later
* Signal Processing Toolbox (for FFT and frequency-domain analysis)
* Statistics and Machine Learning Toolbox (for `tinv` in confidence intervals)

## ⚙️ Getting Started

1. **Clone the repository**

   ```bash
   git clone https://github.com/yourusername/photonic-ode-mrr.git
   cd photonic-ode-mrr
   ```

2. **Open MATLAB** and navigate to the project folder.

3. **Configure user parameters** in `ODE_solver.m` under the **User parameters** section:

   * `order`: system order (1, 2, or 3)
   * `k`: coupling rates per stage (ns⁻¹)
   * `A`: scaling factor (ns → s)
   * Tolerances for coupling (`r_tolerances`), intrinsic loss (`loss_tolerance`), and detuning (`detuning_tolerance`)
   * `N_monte_carlo`: number of Monte Carlo runs
   * Ring geometries (`R`) and effective index (`neff`)

4. **Select input signal** (step or sinusoid) by editing the **Input signal** block.

5. **Run `ODE_solver.m`**

   ```matlab
   ODE_solver
   ```

   * Produces frequency-domain and time-domain plots comparing the ideal ODE solution against the cascaded MRR implementation
   * Executes Monte Carlo simulations and displays envelopes, density plots, and error statistics

## 🔧 Customization

* **Order**: Change `order` to simulate 1st, 2nd, or 3rd-order ODE behavior.
* **Input Signal**: Switch between unit-step (`@(t) C*(t>0)`) and sinusoid (`@(t) C*sin(2*pi*f0*t)`).
* **Tolerances**: Adjust `r_tolerances`, `loss_tolerance`, and `detuning_tolerance` to study robustness.
* **Number of Trials**: Increase `N_monte_carlo` for finer statistical accuracy (at the expense of runtime).

## 🖼️ Results and Visualization

* **Frequency-Domain Plots**: Compare ideal vs. MRR transfer functions; visualize input spectrum and -3 dB bandwidth.
* **Time-Domain Plots**: Show raw input, ODE45 solution, and MRR output, with power comparisons.
* **Monte Carlo Figures**:

  * Percentile envelopes (5–95%) and median of frequency/time responses
  * 2D density heatmap of time-domain outputs
  * RMS error vs. coupling tolerance with confidence intervals

## 📜 License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## 👷🏻‍♂️ Contributing

Feel free to open issues or submit pull requests for improvements, bug fixes, or new features.

## 📫 Contact

For questions or support, contact us at:
* [vittoriopio.cozzoli@mail.polimi.it](mailto:vittoriopio.cozzoli@mail.polimi.it)
* [piervito.creanza@mail.polimi.it](mailto:piervito.creanza@mail.polimi.it)

