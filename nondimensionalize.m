function [M, lam, mu, Rx, tau, X, nu, eta, alfa] ...
    = nondimensionalize(u, v, NEx, m, E, h, nu, t, x, ...
    rho, p, a, gam)
% This function non-dimensionalizes the first 12 input parameters

% rho = 1; p0 = 1/gam; a = 2;
% u0 = 2; 
% NEx = -0.526;
% m = 100;
% E = 72.8E6;
% h = 0.002;
% nu = 0.3;
% endT = 20.0;

D = E*h^3 / (12*(1-nu^2));

M = u * sqrt(rho / (gam*p));
lam = rho*u^2*a^3 / (M*D);
mu = rho*a / m;
Rx = NEx*a^2 / D;
tau = t * sqrt(D / (m*a^4));
X = x/a;
nu; % already nondimensional
eta = h/a;
alfa = atan2(v,u); % atan2(Y,X) <--> tan(alfa) = Y / X
