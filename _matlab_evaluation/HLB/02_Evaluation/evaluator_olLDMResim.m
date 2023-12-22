function evaluator_OL_LDM_resim(segments,config)
%%Function description
%the main task of the function is to generate a .rd5 file that can be used
%in Carmaker out of the provided measurement data

numberOfFiles = size(segments.segments,1);
start = 1;
ending = size(segments.segments(1).segment.q_T0, 1) - 300;
%ending = 30000;
if(numberOfFiles > 1)
    disp("Evaluation is not possible: multiple measurement files.");
else
    segment = segments.segments(1).segment;
    segment.VmcAdas_ldm_egoPoseTheta = segment.VmcAdas_ldm_egoPoseTheta;
    % producing input matrix for the simulation
    input(:,1) = segment.q_T0(start:ending);
    input_corridor.time = segment.q_T0(start:ending);
    input_corridor.signals.values = zeros(4,300,size(segment.q_T0(start:ending), 1));
    input_corridor.signals.dimensions = [4 300];
    input_raw.time = segment.q_T0(start:ending);
    input_raw.signals.values = zeros(4,300,size(segment.q_T0(start:ending), 1));
    input_raw.signals.dimensions = [4 300]; % left border x,y right border x,y 
    for i=start:ending
        corridorXY(i,1:2) = pos_tf2GPS(segment.VmcAdas_ldm_egoPoseX(i),segment.VmcAdas_ldm_egoPoseY(i),segment.VmcAdas_ldm_egoPoseTheta(i),segment.VmcAdas_c0_refline(i));
    end
    corrC1(start:ending+300,1) = atan(pt1_filter(segment.VmcAdas_c1_refline(start:ending+300),1)) + segment.VmcAdas_ldm_egoPoseTheta(start:ending+300);
    corrC2(start:ending+300,1) = movmean(segment.VmcAdas_c2_refline(start:ending+300),100);
    segment.VmcAdas_ldm_egoPoseTheta(1) = segment.VmcAdas_ldm_egoPoseTheta(2);
    for i=start:ending
        distances = ((corridorXY(:,1) - segment.VmcAdas_ldm_egoPoseX(i)).^2 + (corridorXY(:,2) - segment.VmcAdas_ldm_egoPoseY(i)).^2).^0.5;
        startIdx = find(distances == min(distances), 1);
        endIdx = find(distances(startIdx:end) > 150, 1) + startIdx - 1;
        if (endIdx - startIdx > 299)
            endIdx = endIdx - (endIdx - startIdx - 299);
        end
        corridor = corridorXY(startIdx:endIdx,:);

        T = [cos(segment.VmcAdas_ldm_egoPoseTheta(i)) sin(segment.VmcAdas_ldm_egoPoseTheta(i)); -sin(segment.VmcAdas_ldm_egoPoseTheta(i)) cos(segment.VmcAdas_ldm_egoPoseTheta(i))];
        corridorXYEgo = ([corridor(:,1)-segment.VmcAdas_ldm_egoPoseX(i) corridor(:,2) - segment.VmcAdas_ldm_egoPoseY(i)])*T';
        clear corridorXYEgo corrC1 corrC2
        % calculating corridor ego directly from coefficients
        for j=1:151
            corridorXYEgo(j,1) = j-1;
            corridorXYEgo(j,2) = segment.VmcAdas_c0_refline(i) + ...
                segment.VmcAdas_c1_refline(i)*corridorXYEgo(j,1) + ...
                segment.VmcAdas_c2_refline(i)/2*corridorXYEgo(j,1)^2 + ...
                segment.VmcAdas_c3_refline(i)/6*corridorXYEgo(j,1)^3;
            corrC1(j) = segment.VmcAdas_c1_refline(i) + ...
                segment.VmcAdas_c2_refline(i)*corridorXYEgo(j,1) + ...
                segment.VmcAdas_c3_refline(i)/2*corridorXYEgo(j,1).^2;
            corrC2(j) = segment.VmcAdas_c2_refline(i) + ...
                segment.VmcAdas_c3_refline(i)*corridorXYEgo(j,1);
        end
        input_corridor.signals.values(1,1:151,i-start+1) = corridorXYEgo(1:151,1);
        input_corridor.signals.values(2,1:151,i-start+1) = corridorXYEgo(1:151,2);
        input_corridor.signals.values(3,1:151,i-start+1) = corrC1(1:151);
        input_corridor.signals.values(4,1:151,i-start+1) = corrC2(1:151);

        input(i-start+1,2) = 151; %length(corridorXYEgo(1:(endIdx-startIdx+1),1));
    end
    input(:,3) = segment.VmcAdas_ldm_egoPoseX(start:ending);
    input(:,4) = segment.VmcAdas_ldm_egoPoseY(start:ending);
    input(:,5) = segment.VmcAdas_ldm_egoPoseTheta(start:ending);
    tend = input(end,1)-input(1,1);
    tstart = segment.q_T0(start);
    Ts = 0.02;
    assignin('base','input',input);
    assignin('base','input_corridor',input_corridor);
    % for animation
    assignin('base','input_raw',input_raw);
    assignin('base','tend',tend);
    assignin('base','tstart',tstart);
    assignin('base','Ts',Ts);
    
    % starting resim
    rootFolder = config.root(1:end-6);
    cpp2simulink4(rootFolder);
end




end

