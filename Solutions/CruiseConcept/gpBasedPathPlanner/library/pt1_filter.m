function [filtered_y] = pt1_filter(y,c)
filtered_y(1,1) = y(1);
for i=2:length(y)
    filtered_y(i,1) = filtered_y(i-1)*(1-c)+y(i)*c;
end
end

