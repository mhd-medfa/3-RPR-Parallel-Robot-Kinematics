function [pos] = FK_num(q_a, L, z_old)

syms x y L_b phi q real

q11 = q_a(1);
q12 = q_a(2);
q13 = q_a(3);

L_b_num = 0.37;
L_p_num = 0.048;

gamma_vals = [pi/6 5*pi/6 9*pi/6];
% Tleg = Rz(q(1))*Tx(l1)*Tx(q(2))*Tx(l2);


Tbase1 = Rz(pi/6)*Tx(L_b_num);
phi1 = pi/6;
Tbase2 = Rz(5*pi/6)*Tx(L_b_num);
phi2 = 5*pi/6;
Tbase3 = Rz(9*pi/6)*Tx(L_b_num);
phi3 = 9*pi/6;

Ttool1 = Rz(-pi/6)*Tx(-L_p_num);

Ttool2 = Rz(-5*pi/6)*Tx(-L_p_num);

Ttool3 = Rz(-9*pi/6)*Tx(-L_p_num);

gamma_val = gamma_vals(1);
[J1_sym, F1_sym] = FK_sym(q11, L, gamma_val);
% x y L_b gamma phi q real

J1_sym(x ,y, L_b, phi, q) = J1_sym;
F1_sym(x ,y, L_b, phi, q) = F1_sym;

gamma_val = gamma_vals(2);
[J2_sym, F2_sym] = FK_sym(q12, L, gamma_val);

J2_sym(x ,y, L_b, phi, q) = J2_sym;
F2_sym(x ,y, L_b, phi, q) = F2_sym;

gamma_val = gamma_vals(3);
[J3_sym, F3_sym] = FK_sym(q13, L, gamma_val);

J3_sym(x ,y, L_b, phi, q) = J3_sym;
F3_sym(x ,y, L_b, phi, q) = F3_sym;

xd = z_old(1);
yd = z_old(2);
zd = z_old(3);

counter = 0;
z_new = z_old + 5;

while norm(z_new - z_old) > 1e-012
%    disp(norm(z_new - z_old))
   if counter == 0
    z_old = z_old;
   else
    z_old = z_new;
   end
   counter = counter+1;
%    disp(counter)
    %x ,y, L_b, gamma, phi, q
   F1 = F1_sym(xd, yd, L_b_num, phi1, q11); 
   J1 = J1_sym(xd, yd, L_b_num, phi1, q11);
   
   F2 = F2_sym(xd, yd, L_b_num, phi2, q12); 
   J2 = J2_sym(xd, yd,  L_b_num, phi2, q12);
   
   F3 = F3_sym(xd, yd, L_b_num, phi3, q13); 
   J3 = J3_sym(xd, yd, L_b_num, phi3, q13);
   
   JF = double([J1; J2; J3]);
   
   F = double([F1; F2; F3]);
   
   z_new = z_old - pinv(JF)*F;
   
   xd = double(z_new(1));
   yd = double(z_new(2));
   zd = double(z_new(3));
end

pos = [xd, yd, zd];
end

