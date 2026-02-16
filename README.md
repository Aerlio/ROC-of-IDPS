# Robust Optimal Control of Inverter Dominated Power Systems

This repository contains the implementation used to generate the results presented in:

> **T. Ochoa, B. Chaudhuri, M. O’Malley**  
> *Optimal Control for Robust Dynamic Performance in Inverter-Dominated Power Systems*

The source code is written in Julia and MATLAB/Simulink.
This code is shared in the hope that it may be useful to researchers and practitioners working on inverter-dominated power systems. It is provided without warranty and is intended primarily for research and academic use.

## Repository Structure

The entire solution method workflow is orchestrated through `Julia_Main.jl`, which sequentially calls all supporting scripts and manages data flow.

### Julia Scripts

- `Julia_Main.jl` – Executes the full solution method.
- `Julia_1_1_Max.jl` – Maximisation phase ROCP.
- `Julia_1_2_Min.jl` – Minimisation phase ROCP.
- `Julia_1_3_EMT.jl` – Calls MATLAB to perform EMT simulations.
- `Julia_1_4_Plot_series.jl` – Plots time series for a NLP solution and its corresponding EMT-evaluation.
- `Julia_0_0_General_Calls.jl` – General utility functions and MATLAB-Julia connection interface.
- `Julia_0_1_ED.jl` – Economic Dispatch (ED) and Automatic Generation Control (AGC) models.
- `Julia_0_2_Generate_Warmstarts.jl` – Pre-processing step: generates warm-starts for the maximisation phase.
- `Julia_2_1_Plot_results.jl` – Plots various results.

### Matlab Scripts
- `Matlab_1_Generate_Case_Params.m` – Reads test system data.
- `Matlab_2_Ref_Loop.m` – Runs EMT simulation and stores results.
- `Matlab_3_Run_Benchmark_Comparison.m` – Runs EMT simulations over baseline methods using generated disturbance-sequence scenario set.

### Simulink File

- `Simulation.slx` – Contains the DAE-based dynamic model implementation in Simulink.

### Folders

- `CSVs/Case9_3IBR/` - Contains the underlying data for the IEEE 9-bus test system studied in the paper.
- `Generated_Files/`  - Used to store intermediate results throughout the solution method procedure.
 
### Supplementary Material
A detailed mathematical specification of the Robust Optimal Control Problem (ROCP) is provided in the supplementary technical manuscript:

> *Supplementary Mathematical Formulation of the ROCP.pdf*

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
