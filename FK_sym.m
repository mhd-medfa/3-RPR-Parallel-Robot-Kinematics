function [J, F] = FK_sym(q, L, gamma_val)

L1 = L(1);
L2 = L(2);
syms x y L_b phi q real
F = (x - L_b*cos(gamma_val) - L_b*cos(gamma_val+phi))^2 +...
    (y - L_b*sin(gamma_val) - L_b*sin(gamma_val+phi))^2 -...
     (L1 + L2 + q)^2;
% F = (x - L_b*cos(phi) +L_p - L1*cos(phi)*cos(q) )^2 +...
%     (y - L_b*sin(phi)  - L1*sin(phi)*cos(q) )^2 +...
%     (z + L1*sin(q))^2 - (L2^2);
J= [diff(F, x), diff(F, y), 0]; 

end

