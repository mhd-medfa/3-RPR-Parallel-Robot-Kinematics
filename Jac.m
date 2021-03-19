function [J] = Jac(theta, q, phi, L)

L1 = L(1);
L2 = L(2);
L_b_num = 0.37;
gamma_vals = [pi/6 5*pi/6 9*pi/6];

temp1 = - L_b_num*(sin(gamma_vals(1)+phi)+cos(gamma_vals(1)+phi));
temp2 = - L_b_num*(sin(gamma_vals(2)+phi)+cos(gamma_vals(2)+phi));
temp3 = - L_b_num*(sin(gamma_vals(3)+phi)+cos(gamma_vals(3)+phi));

J_z = [ 1 1 temp1;
        1 1 temp2;
        1 1 temp3];

temp4 = (L1+L2+q)*(cos(theta)-sin(theta));
J_theta = diag([temp4, temp4, temp4]);

J =  J_z \ J_theta  ;

end

