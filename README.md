# rush-order-picking
This repository contains the MATLAB script `rush_order_experiment_upload.m` accompanying the article “Real-Time Integration of Rush Orders into Warehouse Order Picking”.

The script generates synthetic problem instances internally; no external input files are required. A fixed random seed (20260926, Twister generator) is used for reproducible instance generation. In the default `paper` mode, the experiment evaluates three scheduling strategies (S0, S1, and S2) across 270 instances.

Requirements: MATLAB with Optimization Toolbox. The MILP formulation is solved using `intlinprog`.

To run the experiment, open and run `rush_order_experiment_upload.m` in MATLAB. Output files are saved automatically in the `rush_order_extended_results` folder next to the script.

Generated outputs include:
- `scenario_inputs.csv`, `order_inputs.csv`, and `picker_inputs.csv`: generated instance data.
- `instance_results.csv`: performance measures and solver information for each instance and strategy.
- `scenario_summary.csv`, `paired_differences.csv`, `paired_summary.csv`, and `deadline_effects.csv`: statistical summaries and comparisons.
- `experiment.mat`: result tables, input tables, and configuration.
- `experiment_settings.txt`: experimental settings.
- `scenario_comparison.png`: comparison plots.

This repository currently provides the code; input and result files are generated when the script is run. Solver runtimes and time-limited MILP results may vary with the MATLAB version and computing environment.
