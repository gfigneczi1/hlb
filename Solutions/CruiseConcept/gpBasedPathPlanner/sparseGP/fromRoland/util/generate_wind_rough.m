function winds = generate_wind_rough(para,type,varargin)

%-------------------------------------------------------------------
%  generate a wind gust at certain time range [t1, t2]
%  
% case-1: constant wind, in line with the wind during traning phase
% case-2: constant-switch wind
% case-3: sine-wind,during [t1,t2]
%-------------------------------------------------------------------
dT = para.dT;
time = para.time;

if isempty(varargin)
    t1 = 0.5*time(end);
    t2 = 1*time(end);
end


const_wind = 1*[2 3 -2];
 amp_wind   =2*[-2 -3 3];

% amp_wind   = -[-3 -4 3];

 switch type
     case 1
         
          winds = ones(length(time),1)*const_wind;
          
     case 2 
         
          winds = [ones(ceil(0.5*length(time)),1)*const_wind;ones(ceil(0.5*length(time))-1,1)*amp_wind];
         
     case 3
         j=1;
%          t1 = 0;

        ws = zeros(length(time),3);

        for i = 0:dT:time(end)
    
            ws(j,:) = const_wind*(0<=i & i<t1)+...
                   (amp_wind.*sin(pi/(t2-t1)*(i-t1))+const_wind)*(i>=t1 & i<=t2)+...
                   const_wind*(i>t2);
             j=j+1;
        end
        winds = ws;            
 end