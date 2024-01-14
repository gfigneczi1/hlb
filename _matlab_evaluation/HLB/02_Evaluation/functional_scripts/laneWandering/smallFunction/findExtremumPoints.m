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
