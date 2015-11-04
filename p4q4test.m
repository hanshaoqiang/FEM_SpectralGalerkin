function e = p4q4test
%P4Q4TEST Test the accuracy of finite element method for project 4
% using an analytical solution.
%
% E = P4Q4TEST
%
% It requires the following MATLAB function signature:
% w = q4(rho0, p0, a, u0, NxE, m, E, h, nu, endT, nstep, paero, nmode, wini, wtini)

nstep   = 10; % 100
N       = 120;
rho0    = 1;
p0      = 1/1.4;
a       = 2;
u0      = 2;
NxE     = -0.526;
m       = 100;
E       = 72.8*10^6;
h       = 0.002;
nu      = 0.3; % by angxiu
endT    = 0.1; % 1
nmode   = 20; % 20
D       = (E*h^3)/(12*(1-nu^2));

x       = repmat(linspace(0,a,N/3+1)',1,nstep);
t       = repmat(linspace(0,endT,nstep),N/3+1,1);

% analytical w(x,t) and derivatives
w       =            sin(pi.*x./a).*sin(pi.*t);
w_xx    = -(pi/a)^2.*sin(pi.*x./a).*sin(pi.*t);
w_xxxx  =  (pi/a)^4.*sin(pi.*x./a).*sin(pi.*t);
w_t     =   pi     .*sin(pi.*x./a).*cos(pi.*t);
w_tt    = - pi^2   .*sin(pi.*x./a).*sin(pi.*t);
Nx      = (E*h*pi^2/(4*a^2)).* (sin(pi.*t).^2) ; % by Shun

paero   = -(D.*w_xxxx - (Nx + NxE).*w_xx + m.*w_tt);
paero   = [zeros(N/3, nstep); paero; zeros(N/3, nstep)]; % by angxiu
wini    = w(:,1);
wtini   = w_t(:,1);

w_numerical     = q4(rho0, p0, a, u0, NxE, m, E, h, nu, endT, nstep, paero, nmode, wini, wtini);
w_analytical    = w;

e_0     = w_analytical - w_numerical;
e       = norm(e_0,'fro')/norm(w_analytical,'fro');
