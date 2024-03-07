function P_matrix = generateParameterSweep(N, maximum, minimum, minStep)
    P = linspace(minimum, maximum, maximum-(N-2)*minStep);
    Nsweep = N*length(P)^2;
    
    for n=1:Nsweep
        for i=1:length(P)
            for j=1:length(P)
                Prem = maximum-P(i)-P(j)-(N-2)*minStep:minStep:maximum-P(i)-P(j);
                P30 = [P(i); P(j); max(10,250 - P1_vector(i)-P2_vector(j)); Prem'];
                P30 = P30/maximum; % Normalization
                P_matrix(:,n) = P30';
            end
        end
    end
end

