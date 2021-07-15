# Thermal runaway avoidance using Hamilton-Jacobi reachability and model predictive control
Andrei Kanavalau, Sanjay Lall

Stanford University

Running the code requires MATLAB and [A Toolbox of Level Set Methods by Ian M. Mitchell](https://www.cs.ubc.ca/~mitchell/ToolboxLS/). Please follow the link for details on how to install and use it.

## Code overview
* System_parameters/System_parameters.m defines the system parameters, which are then used in all other computations/simulations
* Compute_avoid_set/Compute_avoid_set.m computes, saves, and plots the avoid set
* Compute_avoid_set/plot_trajectories/plot_trajectories.m uses the computed avoid set to find points in the vicinity of the boundary. It then computes and plots two trajectories one starting from inside the set, the other --- outside, along with the avoid set.
* MPCs/standard-MPC/standard_MPC.m and MPCs/avoid-MPC/avoid_MPC.m implement simulations of the corresponding MPC schemes
* MPCs/MPC_compare/plot_MPC_results.m plots simulation results
