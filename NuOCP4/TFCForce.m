function [F,U] = TFCForce(q,Pars)
r1 = norm(q-Pars.q1);
r2 = norm(q-Pars.q2);
U1 = -(Pars.gamma/r1);
U2 = -(Pars.gamma/r2);
U = U1 + U2;
F1 = -Pars.gamma*r1^(-3) * (q-Pars.q1);
F2 = -Pars.gamma*r2^(-3) * (q-Pars.q2);
F = F1 + F2;