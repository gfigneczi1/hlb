path = "C:\database\Lexus_steering\csv\sin_70";

csvFiles = dir(fullfile(path, "*.csv"));

steeringRatio = 14.8;
numberOfSpeeds = 3;
numberOfFreq = 4;
numberOfAmp = 2;
freqArray = zeros(numberOfSpeeds,numberOfAmp,numberOfFreq);
amplRatioArray = zeros(numberOfSpeeds,numberOfAmp,numberOfFreq);
phaseShiftArray = zeros(numberOfSpeeds,numberOfAmp,numberOfFreq);

j = 0;

for i=1:size(csvFiles,1)
    data = csvread(fullfile(csvFiles(i).folder, csvFiles(i).name), 1, 0);
    % c1: time, c2: meas, c3: ref, c4: vel
    vx = mean(data(:,4));
    time = data(:,1);
    meas = data(:,2);
    ref = data(:,3);
    dT = mean(diff(time));
    if (contains(convertCharsToStrings(csvFiles(i).name), "step"))
        % step response
    elseif (contains(convertCharsToStrings(csvFiles(i).name), "sine"))
        if (csvFiles(i).name(end-11) == '_')
            freq = str2num(csvFiles(i).name(end-10:end-8))/1000;
        else
            freq = str2num(csvFiles(i).name(end-11:end-8))/1000;
        end
        
        startOfMeas = find(abs(ref) > 0,1) + round(2*1/freq/dT);
                
        %extrema = findExtremumPoints(smooth(ref), []);
        %amplRef = mean(abs(extrema(extrema(:,1)>0,2)));
        amplRef = 0.5*(max(ref(startOfMeas:end))-min(ref(startOfMeas:end)));
        %extrema = findExtremumPoints(smooth(meas), []);
        %amplMeas = mean(abs(extrema(extrema(:,1)>0,2)));
        amplMeas = 0.5*(max(meas(startOfMeas:end))-min(meas(startOfMeas:end)));
                
        zeroCrossingRef = startOfMeas+find(diff(ref(startOfMeas:end)>0)<0);
        zeroCrossingMeas = startOfMeas+find(diff(meas(startOfMeas:end)>0)<0);
        
        phaseShift = mean(zeroCrossingMeas(1:min(length(zeroCrossingMeas),length(zeroCrossingRef))) ...
            -zeroCrossingRef(1:min(length(zeroCrossingMeas),length(zeroCrossingRef))))*dT;
        
        amplRatio = amplMeas/amplRef;
        
        if (vx < 40/3.6)
            j = j+1;
            if (amplRef>30*pi()/180)
                freqArray(1,j,2) = freq;
            else
                freqArray(1,j,1) = freq;
            end
            subplot(2,1,1);
            plot(log(2*pi()*freq), 20*log(amplRatio), 'bo');
            hold on;
            grid on;
            subplot(2,1,2);
            plot(log(2*pi()*freq), -phaseShift*(freq)*2*pi(), 'bo');
            grid on;
            hold on;
        elseif (vx < 60/3.6)
            subplot(2,1,1);
            plot(log(2*pi()*freq), 20*log(amplRatio), 'rx');
            hold on;
            grid on;
            subplot(2,1,2);
            plot(log(2*pi()*freq), -phaseShift*(freq)*2*pi(), 'rx');
            grid on;
            hold on;
        else
            subplot(2,1,1);
            plot(log(2*pi()*freq), 20*log(amplRatio), 'g+');
            hold on;
            grid on;
            subplot(2,1,2);
            plot(log(2*pi()*freq), -phaseShift*(freq)*2*pi(), 'g+');
            grid on;
            hold on;
        end        
    end
    
end