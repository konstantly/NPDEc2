function [F,U] = KeplerForce(q,Pars)
r = norm(q);
U = -(Pars.gamma/r);
F = -Pars.gamma*r^(-3) * q;