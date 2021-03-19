function [q] = IK(Tbase, pos, Ttool, L)

T = eye(4);
T(1:3,4) = pos;

Tleg = Tbase \ T / Ttool;

% check
l1 = L(1);
l2 = L(2);

x = Tleg(1,4);
y = Tleg(2,4);


% Passive angle number 1
q1 = atan2(y , x);

% active angle number 2
r = sqrt(x^2+y^2);
q2 = r - l1 - l2;


q = [q1, q2, 0];

end

