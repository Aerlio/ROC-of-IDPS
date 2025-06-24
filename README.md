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
The entire local reduction method is orchestrated through `1_0_Main.jl`, which sequentially calls all supporting scripts and manages data flow.

### Julia Scripts

- `1_0_Main.jl` – Executes the full local reduction method.
- `1_1_Julia_Max.jl` – Maximisation phase NLP + mesh refinement.
- `1_2_Julia_Min.jl` – Minimisation phase NLP + mesh refinement.
- `1_3_Julia_EMT.jl` – Calls MATLAB to perform EMT simulations.
- `1_4_Julia_Plot_series.jl` – Plots time series for both NLP and EMT results.
- `0_0_Julia_General_calls.jl` – General utility functions and MATLAB-Julia connection interface.
- `0_1_Julia_ED.jl` – Economic Dispatch (ED) and Automatic Generation Control (AGC) models.
- `0_2_Julia_Generate_Warmstarts.jl` – Pre-processing step: generates warm-starts for the maximisation phase.

### Matlab Scripts
- `Matlab_1_Generate_Case_Params.m` – Reads test system data.
- `Matlab_2_Ref_Loop.m` – Runs EMT simulation and stores results.

### Simulink File
- `Simulation.slx` – Contains the IEEE 9 bus test system implementation.

---

## Citation

We kindly request that any publication using this testbed or implementation explicitly acknowledges this repository by citing the aforementioned paper.

---

## Contact

For questions, collaborations, or further information, please contact:

📧 **t.ochoa@ic.ac.uk**

Related publications and presentations are also publicly available at: _[link to be added]_.
