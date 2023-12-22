function data_out = get_velocity_data(file_path, start_of_step, longitudinal_velocity_signal_name)
    AVT_data_table = readtable(file_path);
    long_velocity = table2array(AVT_data_table(:,longitudinal_velocity_signal_name));
    sequence_length = size(long_velocity,1);
    first_second = min([30, sequence_length]);
    velocity_change_difference = max(long_velocity) - min(long_velocity(1:first_second));
    t = linspace(0, sequence_length, sequence_length);
    normalized_long_velocity = (long_velocity- min(long_velocity(1:first_second)))/velocity_change_difference;
    normalized_long_velocity = [zeros(3,1); normalized_long_velocity];
    heaviside_input = heaviside(t-(3+start_of_step));
    heaviside_input = [heaviside_input'; ones(3,1)];
    data_out = iddata(normalized_long_velocity, heaviside_input, 1/30);
    data_out.InterSample = 'foh';
    data_out.Tstart = 0;
end