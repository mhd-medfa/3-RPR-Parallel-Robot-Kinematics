function [T, T1, T2, T3] = FK_test(Tbase,q,Ttool, L)

l1 =L(1);
l2 = L(2);

Tleg = Rz(q(1))*Tx(l1)*Tx(q(2))*Tx(l2);
T1 = Tbase ;
T2 = T1 * Rz(q(1))*Tx(l1) ;
T3 = T2 * Rz(q(1))*Tx(l1)*Tx(q(2))*Tx(l2);
T = Tbase * Tleg * Ttool;

end