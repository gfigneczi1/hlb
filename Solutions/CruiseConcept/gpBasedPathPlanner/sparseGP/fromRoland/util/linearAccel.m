function [ax, ay, az] = linearAccel(RotZYX,Ftotal,para)
   
    
    gz = [0 0 para.g];
    Fz = [0 0 -Ftotal];
    lin_accel = gz' + (1/para.m)*RotZYX*Fz';
    ax = lin_accel(1) ;
    ay = lin_accel(2);
    az = lin_accel(3);
    
end