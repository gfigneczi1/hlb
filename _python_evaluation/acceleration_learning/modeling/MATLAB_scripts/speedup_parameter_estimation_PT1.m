function [K, T] = speedup_parameter_estimation_PT1(speedup_dataset)
    mod = idproc('P1');
    mod.Structure.Kp.Value = 1;
    mod.Structure.Kp.Free = false;
    mod.Structure.Tp1.Minimum = 1;
    mod.Structure.Tp1.Maximum = 10000;
    speedup_model = procest(speedup_dataset, mod);
    K = speedup_model.Kp;
    T = speedup_model.Tp1;
end