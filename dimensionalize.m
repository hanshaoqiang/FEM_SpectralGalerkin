function [u, v, NEx, m, E, h, nu, t, x] ...
    = dimensionalize(M, lam, mu, Rx, tau, X, nu, eta, alfa,...
    rho, p, a, gam)
% This function dimensionalize the first 12 input parameters

u = M / sqrt(rho / (gam*p));
m = rho*a / mu;
x = X*a;
nu;
h = a*eta;
v = u*tan(alfa);

D = rho*u^2*a^3 / (lam*M);

E = D*12*(1-nu^2) / h^3;
NEx = Rx*D / a^2;
t = tau  / sqrt(D / (m*a^4));