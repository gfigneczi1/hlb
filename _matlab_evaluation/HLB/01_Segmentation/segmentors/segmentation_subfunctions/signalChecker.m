function [ mandatoryStatus, testStatus ] = signalChecker( rawData, mandatorySignals, testSignals )
% This function extracts the available fields in the rawData struct.
% Then it checks iff all mandatory signals is available or not.
% If yes, it continues to check if all test data is available.
% If mandatory signals are not available, mandatory status and test status
% both will be set to 0. If mandatory signals are available, their status
% is set to 1, if test data is also available it is set to 1.
mandatoryStatus = 1;
testStatus = 1;

signals = fieldnames(rawData);
for i=1:length(mandatorySignals)
    if (any(signals==mandatorySignals(i))) == 0
        mandatoryStatus = 0;
        testStatus = 0;
        break;
    end
end
if (mandatoryStatus == 1)
    for i=1:length(testSignals)
        if (any(signals==testSignals(i))) == 0
            testStatus = 0;
            break;
        end
    end
end


end

