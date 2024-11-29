function [RotZYX, R_3] = rotBodytoWorld(phi, theta, psi)  % 3-2-1 

Rx_phi  = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
Ry_theta= [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
Rz_psi  = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];

RotZYX = Rz_psi*Ry_theta*Rx_phi;

R13 = RotZYX(1,3);
R23 = RotZYX(2,3);
R33 = RotZYX(3,3);

R_3 = [R13 R23 R33];  % row 1*3

end