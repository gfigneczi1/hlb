function [snippetsPosOffs,snippetsNegOffs, extrema, kpisPosOffs, kpisNegOffs] = functional_laneWanderingSnippetting(offsetError, P0_ma, segment_m, relevantPoints, parameters)

thd = std(offsetError(relevantPoints==1))*0.5;
interactionPoints = abs(offsetError) > thd;
interactionPoints(relevantPoints ~=1) = false;
interactionPoints = morphologyOpen(interactionPoints, 10);
interactionPoints = morphologyClose(interactionPoints, 10);
extrema = findExtremumPoints(offsetError, interactionPoints);

%extrema = calculateExtrema(offsetError,relevantPoints);

%% cut snippets
j = 1;
k = 1;
kpisNegOffs = [];
kpisPosOffs = [];
snippetsNegOffs = [];
snippetsPosOffs = [];

for i=1:size(extrema,1)
    snippetStart = extrema(i,1);
    % snippet start is modified to get the last point, where the thd
    % was exceeded before reaching the extremum
    offsetsBackwards = offsetError(snippetStart:-1:1);
    if (extrema(i,2) > 0)
        % this is a positive offset, searching for first point
        % backwards with negative sign
        thresholdExceed = find(offsetsBackwards <= 0, 1);
    else
        thresholdExceed = find(offsetsBackwards >= 0, 1);
    end

    snippetStop = snippetStart+find(abs(offsetError(snippetStart:end)) < 0.02,1); % the first point where c0 reaches its average after extremum is detected
    if (~isempty(thresholdExceed))
        snippetStart = snippetStart-thresholdExceed;
    end
    if (isempty(snippetStop))
       break;
    end
    compensationAborted = false;
    if (i<size(extrema,1))
        if (snippetStop > extrema(i+1,1))
            % this is the situation, where something happened before
            % getting back to the settle offset, therefore we simply
            % neglect this compensation
            compensationAborted = true;
        end
    end
    if (compensationAborted)
    else
        snippetStop = min(snippetStop, length(offsetError));
        snippetStart = max(snippetStart,1); % adding some init phase
        snippetStopExtended = min(snippetStop + (snippetStop-snippetStart)*2, length(offsetError));
        if (true)
            if (extrema(i,2) > 0)
                % the offset extremum is above the straight line offset
                % of the driver
                snippetsPosOffs(j).name = parameters.id;
                snippetsPosOffs(j).snippet = segment_m(snippetStart:snippetStopExtended, :);
                snippetsPosOffs(j).indexes = parameters.indexes;
                snippetsPosOffs(j).offsetError = offsetError(snippetStart:snippetStopExtended);
                snippetsPosOffs(j).P0_ma = P0_ma(snippetStart:snippetStopExtended);
                snippetsPosOffs(j).relevantPoints.snippetStart = snippetStart;
                snippetsPosOffs(j).relevantPoints.snippetStop = snippetStop;
                snippetsPosOffs(j).relevantPoints.snippetStopExtended = snippetStopExtended;
                snippetsPosOffs(j).relevantPoints.offsetExtremaLocation = extrema(i,1);
                % for positive offset the following quantities are
                % calculated&saved
                % - minimum absolute relative lateral acceleration
                % (right side acceleration)
                % - minimum absolute relative lateral jerk
                % (right side jerk)
                % - minimum absolute curvature change (right side
                % curvature gradient)
                % - maximum lateral offset to the left side
                kpisPosOffs(j, 1:5) = [min(segment_m(snippetStart:snippetStop, parameters.indexes.ayRel)) ...
                    min(segment_m(snippetStart:snippetStop, parameters.indexes.jyRel)) ...
                    min(segment_m(snippetStart:snippetStop, parameters.indexes.vehicleCurvatureChange)) ...
                    min(segment_m(snippetStart:snippetStop, parameters.indexes.TTCL)) ...
                    max(offsetError(snippetStart:snippetStop))];
                 j= j+1;
            else
                snippetsNegOffs(k).name = parameters.id;
                snippetsNegOffs(k).snippet = segment_m(snippetStart:snippetStopExtended, :);
                snippetsNegOffs(k).indexes = parameters.indexes;
                snippetsNegOffs(k).offsetError = offsetError(snippetStart:snippetStopExtended);
                snippetsNegOffs(k).P0_ma = P0_ma(snippetStart:snippetStopExtended);
                snippetsNegOffs(k).relevantPoints.snippetStart = snippetStart;
                snippetsNegOffs(k).relevantPoints.snippetStop = snippetStop;
                snippetsNegOffs(k).relevantPoints.snippetStopExtended = snippetStopExtended;
                snippetsNegOffs(k).relevantPoints.offsetExtremaLocation = extrema(i,1);
                kpisNegOffs(k, 1:5) = [max(segment_m(snippetStart:snippetStop, parameters.indexes.ayRel)) ...
                    max(segment_m(snippetStart:snippetStop, parameters.indexes.jyRel)) ...
                    max(segment_m(snippetStart:snippetStop, parameters.indexes.vehicleCurvatureChange)) ...
                    min(segment_m(snippetStart:snippetStop, parameters.indexes.TTCL)) ...
                    min(offsetError(snippetStart:snippetStop))];
                k = k+1;
            end
        end
    end
end 

end

function extrema = findExtremumPoints(y, inters)
    dy = diff(y);
    dy = [dy; dy(end)];
    posGrad = dy>0;
    negGrad = dy<0;
    maxima = diff(posGrad)<0;
    maxima = [maxima; false];
    minima = diff(negGrad)<0;
    minima = [minima; false];
    extrema = [minima|maxima y];
   % extrema are the global extrema. If inters is empty, this is returned.
   % if not, then inters are considered only.
   if (~isempty(inters))
       % 1 inter means a section of following points
       i = 1; j = 1;
       while i<=length(inters)
           if(inters(i)==1)
               interStop = find(inters(i:end)==0,1);
               interStop = i+interStop-1;
               if (~isempty(interStop))
                   if (mean(y(i:interStop-1)) < 0)
                       [ymin, xmin] = min(y(i:interStop-1));
                       extremaSimplified(j,1:2) = [i+xmin ymin];
                       j = j+1;
                   else
                       [ymax, xmax] = max(y(i:interStop-1));
                       extremaSimplified(j,1:2) = [i+xmax ymax];
                       j = j+1;
                   end
                   i = interStop;
               else
                   break;
               end
           else
               i = i+1;
           end
       end
       extrema = extremaSimplified;
   end
end
