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

## Installing coinhsl library
### Install anaconda
https://docs.anaconda.com/free/anaconda/install/linux/
After that, install the libfortran library
```
conda install -c anaconda libgfortran
```
If that fails (cannot locate the package), run this first:
```
conda config --append channels conda-forge
```
### Download the coinhsl library (license required)
Download the coinhsl-2023.11.17.tar.gz, not the binaries.
After that, extract the tar file and cd to the coinhsl-xxx folder.
Here run the following:
```
sudo apt install libblas-dev liblapack-dev libmetis-dev
```
```
sudo apt install meson
```
```
meson setup builddir --buildtype=release --prefix=/PATH/to_your/coinhsl_folder
```
If you get error for the pkgconfig.relocatable, remove this option from the meson.build file (see line 3, default_options list).
```
meson compile -C builddir
```
```
meson install -C builddir
```
If everything went successfully, follow these steps: https://www.do-mpc.com/en/latest/installation.html#linux (some of these steps are already compoleted)

