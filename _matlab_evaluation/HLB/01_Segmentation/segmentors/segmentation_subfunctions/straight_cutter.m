function [straightIndicator] = straight_cutter(rawData)
% The logic is to have an indicator vector which describes a given scenario
% feature (e.g. straight, curve, highway...etc).
signals = fieldnames(rawData);
N = size(rawData.(signals(1)),1);
straightIndicator = zeros(N,1);

if ~isempty(find(signals=="c2"))
end
end

