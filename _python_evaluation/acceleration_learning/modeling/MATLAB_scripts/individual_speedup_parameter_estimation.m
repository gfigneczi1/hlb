function speedup_parameters = individual_speedup_parameter_estimation(window_directory_path,speedup_index, model_name, time_signal_name,longitudinal_velocity_signal_name)
    file_path = strcat(window_directory_path, "\speed_up", string(speedup_index), ".csv");
    speedup_dataset = get_velocity_data(file_path, 0, longitudinal_velocity_signal_name);
    if model_name == "PT1"
        [K, T] = speedup_parameter_estimation_PT1(speedup_dataset);
        speedup_parameters = [K, T];
    elseif model_name =="CA_PT2"
        speedup_parameters = CA_PT2_optimization(file_path, time_signal_name, longitudinal_velocity_signal_name);
    else
        [K, w0, D] = speedup_parameter_estimation(speedup_dataset);
        speedup_parameters = [K, w0, D];
    end
end