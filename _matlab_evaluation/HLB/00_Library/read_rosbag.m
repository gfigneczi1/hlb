function [topics] = read_rosbag(filePath, savePath)
%READ_ROSBAG Summary of this function goes here
%   Detailed explanation goes here
if (isempty(filePath))
    disp('no measurement file inside');
    topics = [];
else
    bag = rosbag(filePath);
    config_name = fullfile(savePath, 'config.json');
    fid = fopen(config_name);
    config_raw = fread(fid,inf);
    str = char(config_raw');
    fclose(fid);
    config = jsondecode(str);
    topics = cell(length(config.signal_list), 1);
    for i=1:length(config.signal_list)
        topics{i,1} = config.signal_list{i,1}{2,1}{1,1};
    end
    messagesWithoutHeaders = [];
    for i=1:length(topics)
        matFileName = 'topic' + string(i) + '.mat';
        path = fullfile(savePath, matFileName);
        bSel = select(bag,'Topic',topics{i});
        msgCellArrays{i} = readMessages(bSel,'DataFormat','struct');
        if (~isfield(msgCellArrays{i}{1},"Header"))
            disp('A message without header is found');
            messagesWithoutHeaders = [messagesWithoutHeaders i];
        end
    end
    % Checking, if there are any messages without a header
    if (~isempty(messagesWithoutHeaders))
        % searching for proper length header
        for i=1:length(messagesWithoutHeaders)
            properSize = size(msgCellArrays{messagesWithoutHeaders(i)},1);
            for j=1:length(msgCellArrays)
                if (j~=messagesWithoutHeaders(i))
                    if (size(msgCellArrays{j},1) == properSize)
                        disp(strcat(topics(messagesWithoutHeaders(i)),' topic is extended with header artificially.'));
                        for k=1:properSize
                            msgCellArrays{messagesWithoutHeaders(i)}{k}.Header = msgCellArrays{j}{k}.Header;
                        end
                        break;
                    end
                end
            end
        end
    end
    for i=1:length(msgCellArrays)
        if (length(topics{i,1}) == 18)
            if (all(topics{i,1} == '/gps/duro/time_ref'))
                for k=1:size(msgCellArrays{i},1)
                    % time is formed from stamps
                    t(k) = double(msgCellArrays{i}{k}.Header.Stamp.Sec);
                     msgCellArrays{i}{k}.Header.TimeStamp = double(msgCellArrays{i}{k}.Header.Stamp.Sec);
                end
            end
        end
        matFileName = 'topic' + string(i) + '.mat';
        path = fullfile(savePath, matFileName);
        msgCellArray = msgCellArrays{i};
        save(path, 'msgCellArray');
    end
end
end

