function force = indoor_effect_2D(x,y)

% x = xtrans(1);
% y = xtrans(2);
% z = xtrans(3);

%% wall effect

left_wall =  0.1*exp(-0.8*x);
right_wall = 0.1*exp(0.8*x);

front_wall = 0.1*exp(0.8*y);
back_wall = 0.1*exp(-0.8*y);

force = left_wall+ right_wall+ front_wall+ back_wall; % column output


% %% ground & ceiling effect
% 
% ground_effect = exp(-1*z);
% ceiling_effct = exp(1*(z-3))/25;

end