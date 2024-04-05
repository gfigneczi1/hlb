# Model predictive control using do-mpc python library
## Install the requirements
```
pip install -r requirements.txt
```
## Run the simulation
```
python3 nmpc_path_tracking.py
```
## Read the docs
https://www.do-mpc.com/en/latest/

## Parameters
### MPC parameters
```
self.pred_horizon_mpc # Prediction horizon [k]
self.t_step_mpc # Time step of the controller [s]
self.pred_horizon_mpc * self.t_step_mpc # Prediction time [s]
```
### Changing the solver
```
setup_mpc = {
        'n_horizon': self.pred_horizon_mpc,
        'n_robust': 0,
        'open_loop': 0,
        't_step': self.t_step_mpc,
        'state_discretization': 'collocation',
        'store_full_solution': True,
        # Use 'MA27' linear solver in ipopt for faster calculations (current 'mumps'):
        'nlpsol_opts': {'ipopt.linear_solver': 'mumps'} # You have to change the 'mumps' to 'MA27' (not license free)
        }
```
### Cost weights
```
self.alfa
self.beta
self.kappa
```
### Changing the input penalization weight
```
self.mpc.set_rterm(steering_angle=20) # Change the value inside init_controller() function
```
### Simulation parameters
```
self.simulation_steps # Number of iterations in the do_closed_loop() function
self.t_step_simulator
```