# Model Predictive Control at University of Gyor
This solution proposes a self-designed model predictive control (MPC) structure for path follow tasks of real vehicles.
The solution is developed and tested using a real vehicle.

# Problem Formulation
The MPC is requested to accomplish the following maneuvers:
- lane (path) following, where the path is planned by a superordinate component, which is not part of the MPC logic
- desired speed range: 20 kph to 130 kph, with occasional stop&go situation; note: stop&go can be handled by separate control logic
- lateral control (mainly), possible extension with longitudinal control logic
- maximum lateral acceleration of +/- $3.5 m/s^2$

Available infrastructure:
- MATLAB code (simple vehicle models, e.g., kinematic and dynamic bicycle models) together with linear MPC code [1]
- LGSVL simulation, where the four-wheel model of the test Lexus vehicle is available
- ROS framework, where path planning and actuator interfaces are available (ROS2)
- Gazebo simulation framework, where the model of simlpe robots are available + waypoint planner is available
\
<img src="arch.png" alt="architecture" width="800"/>
\
The proposed simple architecture is shown on the above figure. The vehicle control includes the proposed MPC logic.
Low level control consists of two parallel components:
- steering angle control: delivered by the vehicle manufacturer, can be driven by CAN signals,
- speed control: developed by the university research team, which requests the target speed from the vehicle control then translates it acceleration and brake pedal signals (both normalized to range of 0-1). The control is a PID control for both braking and acceleration, including simple anti-windup logic (integrator switch of at saturation) and limitation on both acceleration and jerk for smooth drive-off and stopping. The PID parameters were tuned empirically considering dynamic and comfort requirements. However, the main application area of the developed MPC is lane following, hence speed control is less relevant. Goal is to set a constant maneuver speed and keep it with best accuracy during lane keeping. However, higher lateral acceleration may cause oscillation in the longitudinal speed, in which case the longitudinal control may need further revision.\
\
The univeristy research team uses Autoware Auto packages for various components of the above architecture. Hence, the first solution includes using the MPC proposed by Autoware. Link: https://autowarefoundation.github.io/autoware.universe/main/control/mpc_lateral_controller/
\
This controller is integrated and tuned to the Lexus vehicle, and may be used for reference of comparison in the future. However, our team's goal is to develop an own control logic which is used in various vehicles.

## Integration platforms
The developed algorithm is desired to be used on various platforms:
- Lexus test vehicle
- Formula Student Driverless vehicle (high speed handling)
- F1 1/10 model vehicle (high speed handling of RC vehicle)
- optional: Szenergy Urban Autonomous Challenge vehicles

All the above platforms shall use ROS2 framework.

## MPC interfaces
As mentioned above, path planner component provides the path to follow. This can happen in various representations:
- way point list (array of path points, with predefined even/uneven step sizes, usually given in the vehicle frame or a global frame, each point holds a pose (position+orientation)),
- path parameters (e.g., polynomial coefficients or other curve representations)
\
Additionally the target speed is provided by the speed planner. This may be directly transferred to the low level control or can be further adjusted by the vehicle control. Initially, it is desired to simply feed through the target velocity for simplification.

<img src="mpc.png" alt="mpc" width="600"/>\
The output of the MPC shall be the road wheel angle (or steering angle), which is then delivered to the low level steering angle control. 

## Vehicle Model
Considering the speed and the lateral acceleration requirements, dynamic model is proposed to be used. The vehicle parameters are either available or to be measured. The following parameters may be necessary to properly define the model (list may be incomplete):
- track width and longitudinal axle distance,
- mass,
- steering angle / road wheel angle ratio or look up table,
- possible: dynamic model parameters of the steering system from steering angle to road wheel angle, modelled by a PT2 system -> parameters are planned to be identified in 2024/Q1,
- tire slip characteristics.

# Concepts, ideas
In the beginning phase the following considerations are proposed:
- dynamic vehicle model shall be used,
- vehicle parameters may not be available / may be inaccurate, therefore parameter estimation component may be needed,
- the vehicle model is non-linear, therefore NMPC shall be used,
- non-linear optimization toolbox is to be found, which provides sufficient runtime frequency AND accuracy,
- besides optimization functions the MPC functions may also be taken-over from a tool-box OR developed on our own,
- later, vehicle model may be extended by steering system dynamic model providing better accuracy (to be validated, compared),
- linear implementation, incl. constrained optimization functions is available on MATLAB level [1].

# Results

# References
### Own publications
[1] Ignéczi G., Horvath E., Pup D.: Implementation of a self-developed model predictive control scheme for vehicle parking maneuvers. The first ISTRC Annual Conference, Tel-Aviv, Izrael, 21, June 2021. \
https://www.researchgate.net/publication/354774900_Implementation_of_a_self-developed_model_predictive_control_scheme_for_vehicle_parking_maneuvers

### Literature of control driver models (Cxx)
[C12]   Haobin Jiang, Huan Tian, and Yiding Hua. Model predictive driver model considering the steering characteristics of the skilled drivers. Advances in Mechanical Engineering, 11(3):1–14, 2019. \
[C13]   Alexander Katriniok, Jan P. Maschuw, et al. Optimal vehicle dynamics control for combined longitudinal and lateral autonomous vehicle guidance. In Proceedings of European Control Conference, pages 974–979, Z¨urich, Switzerland, 2013. \
[C14]   J. Kong, M. Pfeiffer, G. Schildbach and F. Borelli, "Kinematic and Dynamic Vehicle Models for Autonomous Driving Control Design," in IEEE Intelligent Vehicles Symposium, Seoul, Korea, 2015.\
[C15]   H. Ye, H. Jiang, S. Ma, B. Tang and L. Wahab, "Linear model predictive control of automatic parking path tracking with soft constraints," International Journal of Advanced Robotic Systems, pp. 1 - 13, May - June 2019. \
[C16]   O. Pauca, C. F. Curuntu and C. Lazar, "Predictive Control for the lateral and longitudinal dynamics in automated vehicles," in 23rd International Conference on System Theory, Control and Computing, Sinaia, 2019. \
[C17]   L. Wang, Model Predictive Control System Design and Implementation using Matlab, London: Springer - Verlag, 2009, pp. 28 - 84, 324 - 359.\
[C18]   Findeisen, Rolf & Allgöwer, Frank. (2002). An introduction to nonlinear model predictive control. 21st Benelux Meeting on Systems and Control. 

### Proposed toolboxes
- https://tinympc.org
- https://gekko.readthedocs.io
- https://docs.acados.org/about_acados/index.html#:~:text=acados%20is%20a%20software%20package,Moritz%20Diehl.
- https://forces.embotech.com/Documentation/_static/FORCESPRO.pdf (for professional application, later phase) 
- https://www.do-mpc.com/en/latest/

### Connecting Repositories
- Path planner using Linear Driver Model (LDM) [6] developed in ROS (obsolete) and ROS2 (maintained)\
  https://github.com/jkk-research/laneKeepSystem
- Szenergy team resources, incl. possible use case of the MPC\
  https://github.com/szenergy/autonomous_master_repo (private)
- Lexus bring up\
  https://github.com/jkk-research/lexus_bringup 