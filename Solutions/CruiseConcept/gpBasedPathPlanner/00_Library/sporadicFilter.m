function [rawData] = sporadicFilter(rawData, cycles)
counter = 0;
currentIndex = 1;
for i=1:length(rawData)
    if (rawData(i)==0)
        if (counter < cycles)
            rawData(currentIndex:i) = 0;
        end
        counter = 0;
        currentIndex = i;
    else
        counter = min(cycles+1,counter + 1);
    end
end
end
    

