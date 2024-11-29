clear; close all;

global time refNorm Frames

path = "C:\database\Lexus_steering\csv\step";

csvFiles = dir(fullfile(path, "*.csv"));

for i=1:size(csvFiles,1)
data = csvread(fullfile(csvFiles(i).folder, csvFiles(i).name), 1, 0);

ref = data(:,2);
time = data(:,1);
% normalizing output
refNorm = (ref-ref(1))/(ref(end)-ref(1));
vx = mean(data(:,4));
A = data(end,3);

D = 0.5; w = 1; tdelay = 3;

Frames = [];

options = optimset('PlotFcns',@optimplotfval);
x = fminsearch(@objectivefcn1,[D w tdelay]);

fprintf("files %d from %d\n", i, size(csvFiles,1));

myVideo = VideoWriter(strcat("Optimization", num2str(i),".avi"));
myVideo.FrameRate = 12;
open(myVideo);
writeVideo(myVideo,Frames);
close(myVideo);

output(i).name = csvFiles(i).name;
output(i).parameters = x(1:2);
output(i).testCase = [vx A];

end

function f = objectivefcn1(x)
global time refNorm Frames
for i=1:length(time)
    y(i) = model(x, time(i));
end

f = sum((y' - refNorm).^2);

idx = find(time>x(3),1);
plot(time(1:idx(1,1)-1),y(1:idx(1,1)-1), 'color', 'k');
hold on;
plot(time(idx(1,1):end),y(idx(1,1):end), 'color', 'm', 'LineWidth',1.5);
plot(time,refNorm);
title(num2str(f));
grid on;

str = strcat("D:", num2str(x(1)));
text(0.15, 0.8, str, 'Color', 'k', 'BackgroundColor', 'w', 'FontSize', 12);
str = strcat("w0:", num2str(x(2)));
text(0.15, 0.7, str, 'Color', 'k', 'BackgroundColor', 'w', 'FontSize', 12);
str = strcat("tdelay:", num2str(x(3)));
text(0.15, 0.6, str, 'Color', 'k', 'BackgroundColor', 'w', 'FontSize', 12);
xlabel('time(s)');ylabel('speed normalized to settle value');
hold off;

Frames = [Frames getframe(gcf)]; 

end




function Y = model(P, U)
% returns output value based on time = U
D = P(1);
T = 1/P(2);
tc = 0;
alfa = 0;
tdelay = P(3);

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