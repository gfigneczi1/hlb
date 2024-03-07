function f = functional_optimizeNodePointDistances(P3)

global P3in opt options
P3in = P3;

[traj, ref, ~, ~, ~, ~, ~, GT_U, GT_Y]  = functional_driverModelWithPlanner(options.parameters.P0, 0);

f = evaluation_lossCalculationTrajectoryPlanner(traj(:,:),ref(:,:),GT_U, GT_Y, opt, P3in);

end

