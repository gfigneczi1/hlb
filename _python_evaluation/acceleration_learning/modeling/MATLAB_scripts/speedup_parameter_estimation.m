function  [K, w0, D] = speedup_parameter_estimation(speedup_dataset)
    mod = idproc('P2U');
    mod.Structure.Kp.Value = 1;
    mod.Structure.Kp.Free = false;
    mod.Structure.Tw.Value = 2;
    mod.Structure.Tw.Minimum = 1;
    mod.Structure.Tw.Maximum = 10000;
    mod.Structure.Zeta.Value = 0.5;
    mod.Structure.Zeta.Maximum = 0.9999;
    speedup_model = procest(speedup_dataset, mod);
    K = speedup_model.Kp;
    w0 = 1/speedup_model.Tw;
    D = speedup_model.Zeta;
end