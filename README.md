# Robust Optimal Control of Inverter Dominated Power Systems

This repository contains the MATLAB/Simulink model and Julia/MATLAB code used to generate the results presented in the paper:

> **T. Ochoa, B. Chaudhuri, M. O’Malley**  
> *Optimal Control for Robust Dynamic Performance in Inverter-Dominated Power Systems (Part I & II)*

The source code is provided in the hope that it will be useful to researchers and practitioners, but it comes without any warranty.

## Repository Structure

- `CSVs/Case9_3IBR/`  
  Contains the underlying data for the IEEE 9-bus test system studied in the paper.

- `Generated_Files/`  
  Used to store results produced by the local reduction method.
  
## Scripts
The entire local reduction method is orchestrated through `0_Main.jl`, which sequentially calls all supporting scripts and manages data flow.

### Julia Scripts

- `0_Main.jl` – Executes the full local reduction method.
- `Julia_1_1_Max.jl` – Maximisation phase NLP + mesh refinement.
- `Julia_1_2_Min.jl` – Minimisation phase NLP + mesh refinement.
- `Julia_1_3_EMT.jl` – Calls MATLAB to perform EMT simulations.
- `Julia_1_4_Plot_series.jl` – Plots time series for both NLP and EMT results.
- `Julia_1_5_Plot_score.jl` – Plots main results of local reduction method.
- `Julia_0_0_General_calls.jl` – General utility functions and MATLAB-Julia connection interface.
- `Julia_0_1_ED.jl` – Economic Dispatch (ED) and Automatic Generation Control (AGC) models.
- `Julia_0_2_Generate_Warmstarts.jl` – Pre-processing step: generates warm-starts for the maximisation phase.

### Matlab Scripts
- `Matlab_1_Generate_Case_Params.m` – Reads test system data.
- `Matlab_2_Ref_Loop.m` – Runs EMT simulation and stores results.

### Simulink File
- `Simulation.slx` – Contains the IEEE 9 bus test system implementation.

---

## Solver Requirements

This project uses the [KNITRO](https://www.artelys.com/solvers/knitro/) nonlinear optimization solver for solving large-scale NLPs.  
> ⚠️ **KNITRO requires a valid license.**  
Ensure you have a valid KNITRO license installed and configured before running the Julia scripts that invoke the solver.

The connection to KNITRO is handled through the Julia interface to `KNITRO.jl`.

---

## Citation

We kindly request that any publication using this testbed or implementation explicitly acknowledges this repository by citing the aforementioned paper.

---
## Contact

For questions, collaborations, or further information, please contact:

📧 **t.ochoa@ic.ac.uk**

Related publications and presentations are also publicly available at: _[link to be added]_.
