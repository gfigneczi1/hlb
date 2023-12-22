function Conv2TimeseriesPostProcess(InputFolder)

matFiles = {};
temp = dir(fullfile(InputFolder,'*.mat'));
for i = 1:length(temp)
    matFiles{end+1} = fullfile(InputFolder, temp(i).name);
end
emptyCell = {};
for i = 1:length(matFiles)
    data = load(matFiles{i});
    if isfield(data, 'q_T0')
        data = parseStruct(data, emptyCell, data, data.q_T0);
        data = rmfield(data, 'q_T0');
        save(matFiles{i}, '-struct', 'data')
    end
end
end

function output = parseStruct(input, list, output, q_T0)


if isnumeric(input)
    temp = getfield(output,list{:});
    tempTimeSeries = timeseries(temp', q_T0');
    output = setfield(output,list{:}, tempTimeSeries);
else
    temp = fieldnames(input);
    for i =1:length(temp)
        list{end+1} = temp{i};
        output = parseStruct(input.(temp{i}), list, output, q_T0);
        list(end)=[];
    end
end


end