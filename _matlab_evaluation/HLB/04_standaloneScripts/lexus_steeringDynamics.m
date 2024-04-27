% this is a simple simulation for inverse dynamics
clear all;
Td = 0.02;
t = 0:Td:5;
x(:,1) = [0;0];
x_ref = 1;
ie = 0;
umax = 1;
e_k1 = 0;
% simulation between 0-5secs
for i=2:length(t)
    e = x_ref-x(1,i-1);
    ie = ie+e*Td;
    de = (e-e_k1)/Td;
    e_k1 = e;
    u = controller(e, ie, de);
    u = max(min(u,umax),-umax);
    u = inverseModel(x(:,i-1),u);
    x(:,i) = model(x(:,i-1),u);
end
% end of simulations
plot(t,x(1,:));

function x = model(x_k1, u)
    T = 0.6; kszi = 1.2; Td = 0.02; K =1;
    A0 = 1+2*kszi*T/Td+T^2/Td^2;
    x(1) = K/A0*(u+2*T^2/Td^2*x_k1(1)-T^2/Td^2*x_k1(2)+2*kszi*T/Td*x_k1(1));
    x(2) = x_k1(1);
end

function y =inverseModel(x_k1,u)
    T = 0.2; kszi = 0.6; Td = 0.02;
    y = x_k1(2)*T^2/Td^2+x_k1(1)*(-2*kszi*T/Td-2*T^2/Td^2)+u*(T^2/Td^2+2*kszi*T/Td+1);
end

function y = controller(e, ie, de)
    % simple PID controller, only for tests
    P = 2;
    I = 0.3;
    D = 0.7;
    y = P*e+I*ie+D*de;
end