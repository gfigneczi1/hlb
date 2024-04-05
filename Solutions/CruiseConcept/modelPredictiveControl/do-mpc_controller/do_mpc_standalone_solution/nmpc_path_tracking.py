import do_mpc
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.ndimage import uniform_filter1d
import math
from casadi import *
from timeit import default_timer
import matplotlib.pyplot as plt

class MPCController:
    def __init__(self) -> None:
        df = pd.read_excel('31_south.xlsx')
        self.x_coordinates = df.iloc[:, 0].values
        self.y_coordinates = df.iloc[:, 1].values
        self.simulation_steps = 11000 # simulation time: simulation_steps * t_step_simulator
        self.ref_road_points_global = np.column_stack((self.x_coordinates, self.y_coordinates))
        model_type = 'continuous'
        self.model = do_mpc.model.Model(model_type)
        # Vehicle Model Parameters 
        self.vx = 18 # [m/s]
        self.Calfa_f = 100000 # [N/rad]
        self.Calfa_r = 100000 # [N/rad]
        self.Lr = 1.8 # [m]
        self.Lf = 1.2 # [m]
        self.m = 1500 # [kg]
        self.Iz = 2000 # [kg*m^2]
        # cost weights 
        self.alfa = np.array([[10, 0], [0, 1]])
        self.beta = 1.0
        self.kappa = 0
        # control parameters
        self.pred_horizon_mpc = 64 #50
        self.t_step_mpc = 0.02 #0.03
        # for dummy testing
        self.y_ref =  np.array([[3], [0]])
        # for road testing (just init, calculate later)
        self.y_ref_vec_global = np.zeros(self.pred_horizon_mpc+4)
        self.orientation_ref_vec_global = np.zeros(self.pred_horizon_mpc+4)
        # set up the model
        self.init_model()
        # set up the controller
        self.mpc = do_mpc.controller.MPC(self.model)
        self.init_controller()
        # set up estimator, simulator
        self.t_step_simulator = 0.005
        self.estimator = do_mpc.estimator.StateFeedback(self.model)
        self.simulator = do_mpc.simulator.Simulator(self.model)
        self.init_simulator()
        # set up initial states
        #self.x0 = np.array([0, 0, self.x_coordinates[0], self.y_coordinates[0], 0])
        self.x0 = np.array([0, 0, 0, 0, 0])

        self.simulator.x0 = self.x0
        self.mpc.x0 = self.x0
        self.estimator.x0 = self.x0
        self.mpc.set_initial_guess()

        # for debug
        self.error_vec = np.zeros((1, 2))

    def init_model(self):
        """
        This function sets up the corresponding ODE to the model
        """
        self.x1 = self.model.set_variable('_x',  'x1', shape=(1,1))
        self.x2 = self.model.set_variable('_x',  'x2', shape=(1,1))
        self.x3 = self.model.set_variable('_x',  'x3', shape=(1,1))
        self.x4 = self.model.set_variable('_x',  'x4', shape=(1,1))
        self.x5 = self.model.set_variable('_x',  'x5', shape=(1,1))

        self.y_ref0 = self.model.set_variable('_tvp', 'y_ref0', shape=(1,1))
        self.theta_ref0 = self.model.set_variable('_tvp', 'theta_ref0', shape=(1,1))

        self.u = self.model.set_variable('_u',  'steering_angle')

        # Algebraic equations
        A = - (self.Calfa_f * cos(self.u) + self.Calfa_r) / (self.m * self.vx)
        B = (-self.Lf * self.Calfa_f * cos(self.u) + self.Lr * self.Calfa_r) / (self.m * self.vx) - self.vx
        C = (-self.Lf * self.Calfa_f * cos(self.u) + self.Lr * self.Calfa_r) / (self.Iz * self.vx)
        D = (-self.Lf**2 * self.Calfa_f * cos(self.u) + self.Lr**2 * self.Calfa_r) / (self.Iz * self.vx)
        E = (self.Calfa_f * cos(self.u)) / self.m
        F = (self.Lf * self.Calfa_f * cos(self.u)) / self.Iz

        # RIGHT HAND SIDE EQUATIONS
        self.model.set_rhs('x1', A * self.x1 + C * self.x2 + E * self.u)
        self.model.set_rhs('x2', B * self.x1 + D * self.x2 + F * self.u)
        self.model.set_rhs('x3', self.vx * cos(self.x5) - self.x1 * sin(self.x5))
        self.model.set_rhs('x4', self.vx * sin(self.x5) + self.x1 * cos(self.x5))
        self.model.set_rhs('x5', self.x2)

        # not used for path tracking
        self.model.set_expression(expr_name='cost1', expr=self.beta * (np.linalg.norm(np.dot((np.array([[self.x4], [self.x5]]) - self.y_ref).T, (np.array([[self.x4], [self.x5]]) - self.y_ref))))) # terminal cost  
        self.model.set_expression(expr_name='cost2', expr=(np.dot((np.array([[self.x4], [self.x5]]) - self.y_ref).T, np.dot(self.alfa, (np.array([[self.x4], [self.x5]]) - self.y_ref))) + self.kappa * self.u ** 2))# stage cost 

        self.model.setup()

    def init_controller(self):
        """
        This function sets up the controller. The cost functions are also defined here.
        """
        setup_mpc = {
        'n_horizon': self.pred_horizon_mpc,
        'n_robust': 0,
        'open_loop': 0,
        't_step': self.t_step_mpc,
        'state_discretization': 'collocation',
        'store_full_solution': True,
        # Use 'MA27' linear solver in ipopt for faster calculations (current 'mumps'):
        'nlpsol_opts': {'ipopt.linear_solver': 'MA57'}
        }
        self.mpc.set_param(**setup_mpc)

        # mterm = model.aux['cost1'] # terminal cost
        # lterm = model.aux['cost2'] # stage cost
        mterm = self.beta * (np.linalg.norm(np.dot((np.array([[self.x4], [self.x5]]) - np.array([[self.y_ref0], [self.theta_ref0]])).T, (np.array([[self.x4], [self.x5]]) - np.array([[self.y_ref0], [self.theta_ref0]]))))) # terminal cost
        lterm = np.dot((np.array([[self.x4], [self.x5]]) - np.array([[self.y_ref0], [self.theta_ref0]])).T, np.dot(self.alfa, (np.array([[self.x4], [self.x5]]) - np.array([[self.y_ref0], [self.theta_ref0]])))) + self.kappa * self.u ** 2 # stage cost

        self.mpc.set_objective(mterm=mterm, lterm=lterm)
        self.mpc.set_rterm(steering_angle=20)
        self.mpc.settings.supress_ipopt_output()

        tvp_template = self.mpc.get_tvp_template()

        def tvp_fun(t_now):
            """
            Setting the time varying parameters based on the current control step. 
            """
            for k in range(self.pred_horizon_mpc + 1):
                if k >= self.pred_horizon_mpc-2:
                    tvp_template['_tvp', k, 'y_ref0'] = self.y_ref_vec_global[-1]
                    tvp_template['_tvp', k, 'theta_ref0'] = self.orientation_ref_vec_global[-1]
                else:
                    y_ref0_interp_fun = interpolate.interp1d(np.array([k*self.t_step_mpc, (k+1)*self.t_step_mpc]), np.array([self.y_ref_vec_global[k+1], self.y_ref_vec_global[k+2]]))
                    theta_ref0_interp_fun = interpolate.interp1d(np.array([k*self.t_step_mpc, (k+1)*self.t_step_mpc]), np.array([self.orientation_ref_vec_global[k+1], self.orientation_ref_vec_global[k+2]]))
                    tvp_template['_tvp', k, 'y_ref0'] = y_ref0_interp_fun(k*self.t_step_mpc)
                    tvp_template['_tvp', k, 'theta_ref0'] = theta_ref0_interp_fun(k*self.t_step_mpc)

            return tvp_template

        self.mpc.set_tvp_fun(tvp_fun)

        # self.mpc.bounds['lower', '_u', 'steering_angle'] = -2
        # self.mpc.bounds['upper', '_u', 'steering_angle'] = 2
        self.mpc.setup()

    def init_simulator(self):
        """
        This function sets up the simulator. The time varying parameters have to be set in a different way
        read more about it in the do-mpc docs.
        """
        self.simulator.set_param(t_step = self.t_step_simulator)
        
        tvp_template = self.simulator.get_tvp_template()
        
        def tvp_fun(t_now):
            tvp_template['y_ref0'] = self.y_ref_vec_global[0]
            tvp_template['theta_ref0'] = self.orientation_ref_vec_global[0]

            return tvp_template

        self.simulator.set_tvp_fun(tvp_fun)

        self.simulator.setup()

    def visualize_results(self):
        '''
        This function visualizes the system states, manipulated variable, and costs over the simulation time,
        also saves the error distance between the vehicle and the reference path, and the reference path vs ego path
        '''
        from matplotlib import rcParams
        plt.figure()
        plt.plot(self.error_vec[:,1], self.error_vec[:,0])
        plt.savefig('error_debug.png')
        plt.figure(2)
        plt.plot(self.x_coordinates[:len(self.mpc.data['_x', 'x4'])], self.y_coordinates[:len(self.mpc.data['_x', 'x4'])])
        plt.plot(self.mpc.data['_x', 'x3'], self.mpc.data['_x', 'x4'], color='r')
        plt.savefig('road_pos_debug.png')
        rcParams['axes.grid'] = True
        rcParams['font.size'] = 18
        fig, ax, graphics = do_mpc.graphics.default_plot(self.mpc.data, figsize=(16,9))
        graphics.plot_results()
        graphics.reset_axes()
        plt.savefig('simu_results.png')
        plt.show()

    def monotonic_safe_interpolation(self, ref_road_points_local, xend, end_point):
        '''
        This function looks for segments in the reference road, where the X coordinate is monotonic,
        currently only uses the first segment, when the road segment is divided into multiple segments.
        '''
        start_indices = [0]

        # Step 2
        increasing = ref_road_points_local[1, 0] > ref_road_points_local[0, 0]
        for i in range(2, len(ref_road_points_local[:end_point])):
            if (ref_road_points_local[i, 0] > ref_road_points_local[i-1, 0]) != increasing:
                # Step 3
                start_indices.append(i-1)
                increasing = not increasing

        # Step 4
        start_indices.append(len(ref_road_points_local[:end_point]))

        # Step 5
        y_ref_vec = np.empty(self.pred_horizon_mpc+1)
        for i in range(len(start_indices) - 1):
            start, end = start_indices[i], start_indices[i+1]
            try:
                interp_fn = interpolate.CubicSpline(ref_road_points_local[start:end, 0], ref_road_points_local[start:end, 1])
            except:
                plt.figure()
                plt.plot(self.error_vec[:,1], self.error_vec[:,0])
                plt.savefig('error_debug.png')
                plt.figure(2)
                plt.plot(self.x_coordinates[:len(self.mpc.data['_x', 'x4'])], self.y_coordinates[:len(self.mpc.data['_x', 'x4'])])
                plt.plot(self.mpc.data['_x', 'x3'], self.mpc.data['_x', 'x4'], color='r')
                plt.savefig('road_pos_debug.png')
            x_ref_vec = np.linspace(0, xend, self.pred_horizon_mpc+1)
            y_ref_vec[start:self.pred_horizon_mpc+1] = interp_fn(x_ref_vec)
            break

        return y_ref_vec

    def calculate_reference_vectors(self):
        """
        This function calculates the reference vector Y, and reference vector Theta in the given prediction horizon.
        The cost functions will use these values later.
        """
        nearest_index = np.argmin(np.sum((self.ref_road_points_global - np.array([self.x0[2], self.x0[3]]).T)**2, axis=1)**0.5)
        dx = self.vx * self.t_step_mpc
        xend = dx * self.pred_horizon_mpc

        rot_matrix = np.array([[cos(self.x0[4]), sin(self.x0[4])], [-sin(self.x0[4]), cos(self.x0[4])]])
        transl_matrix = np.array([self.x0[2], self.x0[3]])
        ref_road_points_local = np.dot((self.ref_road_points_global - transl_matrix.T), rot_matrix.T)

        local_distances = np.sum(ref_road_points_local ** 2, axis=1)**0.5
        end_point = nearest_index + np.argmax(local_distances[nearest_index:] >= xend)
        x_ref_vec = np.linspace(0, xend, self.pred_horizon_mpc+1)
        #interp_fn = interpolate.CubicSpline(ref_road_points_local[:end_point, 0], ref_road_points_local[:end_point, 1])
        #y_ref_vec = interp_fn(x_ref_vec)
        y_ref_vec = self.monotonic_safe_interpolation(ref_road_points_local, xend, end_point)
        orientation_ref_vec = np.diff(y_ref_vec) / np.diff(x_ref_vec)
        orientation_ref_vec = uniform_filter1d(orientation_ref_vec, size=10) # not quite the same result as matlab movmean
        orientation_ref_vec = np.append(orientation_ref_vec, orientation_ref_vec[-1])
        points_ref_global = np.dot(np.column_stack((x_ref_vec, y_ref_vec)), rot_matrix) + transl_matrix.T
        self.x_ref_vec_global = points_ref_global[:,0]
        self.y_ref_vec_global = points_ref_global[:,1]
        self.orientation_ref_vec_global = orientation_ref_vec + self.x0[4]

    def do_closed_loop(self):
        """
        This function is responsible for the closed loop simulation.
        """
        for k in range(self.simulation_steps):
            t1 = default_timer()
            # First calculate the reference vectors
            self.calculate_reference_vectors()
            print(self.x0)
            # Then calculate the manipulated variable (control var)
            u0 = self.mpc.make_step(self.x0)
            # Calculate the error distance for debug purposes
            current_error = np.array([[(float(self.mpc.data['_x', 'x3'][k] - self.x_ref_vec_global[0])**2 + float(self.mpc.data['_x', 'x4'][k] - self.y_ref_vec_global[0])**2)**0.5, k]])
            self.error_vec = np.concatenate((self.error_vec, current_error), axis=0)
            # Apply the manipulated variable to the system
            y_next = self.simulator.make_step(u0)
            # print(y_next)
            # Measure the system states after the control cycle
            self.x0 = self.estimator.make_step(y_next)
            t2 = default_timer()
            print(f'One control loop took {t2 - t1} time.')
   
if __name__ == '__main__':
    mpc_instance = MPCController()
    mpc_instance.do_closed_loop()
    mpc_instance.visualize_results()