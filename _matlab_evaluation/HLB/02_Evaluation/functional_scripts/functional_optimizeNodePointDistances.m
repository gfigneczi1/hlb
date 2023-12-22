function f = functional_optimizeNodePointDistances(P3)

global P3in opt
P3in = P3;
P0 = zeros(21,1); % dummy model parametrization

[traj, ref, ~, ~, ~, ~, ~, GT_U, GT_Y]  = functional_driverModelWithPlanner(P0, 0);

f = evaluation_lossCalculationTrajectoryPlanner(traj(:,:),ref(:,:),GT_U, GT_Y, opt, P3in);

end

