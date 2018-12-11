# Coded Caching Simulator for MISO Broadcast Channels.

This simulator can be used to generate and verify, e.g., the numerical results provided in article [Multi-antenna Interference Management for Coded Caching](https://arxiv.org/abs/1711.03364).

The simulator has been written in MATLAB and uses [CVX](http://www.cvxr.com) to solve convex optimization problems.

## Usage

The simulator relies on [CVX](http://www.cvxr.com) to solve convex problems. Either the cvx folder needs to be placed in the simulator root or CVX needs to be globally configured.

The simulation scenarios are stored in the sims folder under individual subfolders.

### Example: 6 user scenario with varying alpha parameter

`cd sims/K6-alpha`

Run for all different precoder designs

`simulate_subset_2.m`

`simulate_subset_3.m`

`simulate_subset_4.m`

`simulate_subset_5.m`

`simulate_subset_6.m`

Plot the results

`plotsim`
