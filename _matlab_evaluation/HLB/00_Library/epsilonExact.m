function [epsExact,epsExactIndx]=epsilonExact(pointsAll, points)
%Find the largest perpendicular distance from a simplified curve (points) to
%the original curve (pointsAll).
%Return also the index of the point with the largest perpendicular distance.
%INPUT:
%      pointsAll: points of the original curve, n x [x, y, id]
%      points:    points of the simplified curve, m x [x, y, id], m <= n
%OUTPUT:
%       epsExact: perpendicular distance of the point with the largest distance to the simplified curve.
%       epsExactIndx: index of the point pointsAll(epsExactIndx) with the largest perpendicular distance.
%EXAMPLE: 
%         x=[2;3;4]; y=[0;1;0]; id=[1;2;3]; pointsAll=[x,y,id];
%         x=[2;4]; y=[0;0]; id=[1;3]; points=[x,y,id];
%         [epsExact,epsExactIndx]=epsilonExact(pointsAll, points)
numpointsOut=size(points,1);% number points of simplified curve
epsExact=0;
for k=1:numpointsOut-1
  id1=points(k,3);%identifier of 1st point of line segment
  id9=points(k+1,3);%identifier of last point of line segment
  indx1=find(pointsAll(:,3)==id1);%index of pointsAll with id1
  pointsCur=pointsAll(indx1:find(pointsAll(:,3)==id9),:);% all points of the line segment
  if size(pointsCur,1)<3
    %Adjacent pointsCur
    continue
  end
  
  % Perpendicular distances d of the pointsCur Pp to the line segment
  % Pp are all pointsCur within line segment range
  Pp = pointsCur(2:end-1,1:2);
  x1 = pointsCur(1,1);
  y1 = pointsCur(1,2);
  dx = pointsCur(end,1)-x1;% x2-x1
  dy = pointsCur(end,2)-y1;% y2-y1
  den = sqrt(dx^2+dy^2);% denominator
  if den > eps
    % segment is a line
    d = abs(dx*(y1-Pp(:,2))-(x1-Pp(:,1))*dy)/den;
  else
    % segment is a point (line with zero length)
    d = sqrt((x1-Pp(:,1)).^2+ (y1-Pp(:,2)).^2);
  end
  [pDistsMax,pDistsMaxIdx] = max(d);%[largest perpend. dist. of this line segment, index of it]
  
  if pDistsMax>=epsExact
    epsExact=pDistsMax;%largest epsilon so far
    epsExactIndx=pointsAll(indx1+pDistsMaxIdx,3);
  end
end