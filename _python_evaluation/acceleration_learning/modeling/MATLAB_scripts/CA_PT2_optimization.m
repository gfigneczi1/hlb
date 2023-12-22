function parameters_CA_PT2 = CA_PT2_optimization(file_path, time_signal_name, longitudinal_velocity_signal_name)
    global time refNorm
    close all;
    AVT_data_table = readtable(file_path);
    ref = table2array(AVT_data_table(:,longitudinal_velocity_signal_name));
    time = table2array(AVT_data_table(:,time_signal_name)) - table2array(AVT_data_table(1,time_signal_name));
    % normalizing output
    refNorm = (ref-ref(1))/(ref(end)-ref(1));
    
    D = 0.5; w = 1; tc = 3; alfa = 0.125; tdelay = 0.5;
    
    options = optimset('PlotFcns',@optimplotfval);
    x = fminsearch(@objectivefcn1,[D w tc alfa tdelay], options);
    
    for i=1:length(time)
        y(i) = model(x, time(i));
    end
    parameters_CA_PT2 = x;
    close all;
end

function f = objectivefcn1(x)
global time refNorm
for i=1:length(time)
    y(i) = model(x, time(i));
end

f = sum((y' - refNorm).^2);

%plot(time,y);
%hold on;
%plot(time,refNorm);
%title(num2str(f));
%grid on;
%hold off;

end

function Y = model(P, U)
% returns output value based on time = U
D = P(1);
T = 1/P(2);
tc = P(3);
alfa = P(4);
tdelay = P(5);

v0 = alfa*tc; G = 1;

if (U<tdelay)
    Y = 0;
elseif (U<(tc+tdelay))
    Y = alfa*(U-tdelay);
else
    t = U-tc-tdelay;
    if (D>1)
        l1 = (-D+(D^2-1)^0.5)/T;
        l2 = (-D-(D^2-1)^0.5)/T;
        c1 = (v0-G)*(l2-alfa)/(l2-l1);
        c2 = (alfa-l1*(v0-G))/(l2-l1);
        Y = c1*exp(l1*t) + c2*exp(l2*t) + G;
        %Y = l2/(l1-l2)*exp(l1*t) + l1/(l2-l1)*exp(l2*t) + 1;
    elseif (D == 1)
        l = -1/T;
        c1 = v0-G;
        c2 = alfa - l*(v0-G);
        Y = c1*exp(l*t) + c2*t*exp(l*t) + G;
        %Y = 1-exp(-t/T)*(1+t/T);
    else
        A = -D/T;
        B = ((1-D^2)^0.5)/T;
        c1 = v0-G;
        c2 = (alfa-A*(v0-G))/B;
        Y = exp(A*t)*(c1*cos(B*t) + c2*sin(B*t)) + G;
        %Y = 1-exp(-D*t/T)*(cos(B*t)+D/(T*B)*sin(B*t));
    end
    %Y = K*Y + alfa*tc;
end
end