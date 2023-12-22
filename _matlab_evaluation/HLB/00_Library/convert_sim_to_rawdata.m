function [meas_out] = convert_sim_to_rawdata(sim_in)
%CONVERT_SIM_TO_RAWDATA Summary of this function goes here
%   convert model_in from VMC sim to a mat file for the matlab eval
input = load(sim_in);
fields = fieldnames(input.model_in);
len = size(input.model_in.X_abs.Data,1);

for i=1:length(fields)
    field = fields{i};
    meas_out.(field) = (input.model_in.(field).Data)';
    if length(meas_out.(field)) == 1
        meas_out.(field) = ones(1,len) .* meas_out.(field);
    end
end
save('validator_62_20221219_.mat','-struct','meas_out');
end

