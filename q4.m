function [w,varargout] = q4(rho0, p0, a, u0, NEx, m, E, h, nu,...
    endT, nstep, paero, nmode, wini, wtini)

%% parameter
gam = 1.4;
v0 = 0;

D = E*h^3 / (12*(1-nu^2));
dtaudt = sqrt( D / (m*a^4) ); % dtau/dt

%% derived parameters
N = size(paero,1) - 1; % number of horizontal (x direction) elements across the flow field
x = linspace(0,a,N/3+1)'; % coordinates of nodes on the plate
dx = 3*a/N; % element size
% CFL = u0*endT / (nstep * dx)

dt = endT / (nstep - 1);
dofs = 2*nmode; % global degrees of freedom

%% parse initial condition
% % convert from wini & wtini to Aini and Bini by evaluating inner product of
% % wini or wtini and basis functions.
% Phi_avg = 0.5*( Phi(1:end-1,:) + Phi(2:end,:) );
% 
% wini_avg = 0.5*( wini(1:end-1) + wini(2:end) );
% wtini_avg = 0.5*( wtini(1:end-1) + wtini(2:end) );

Phi = basis(nmode,x,a); % basis functions evaluated at x
Npt = size(wini,1); % number of nodes on plate; used for discrete Fourier transform
Aini = Phi' * wini * 2 / (h*(Npt-1));
Bini = Phi' * wtini * 2 / (h*(Npt-1)*dtaudt);

U0 = zeros(dofs,1);
U0([1:2:end-1]) = Aini;
U0([2:2:end]) = Bini;

%% paramters wrapup
param = struct;
param.gam = gam;
param.rho = rho0;
param.p0 = p0;
param.a = a;
param.u0 = u0;
param.v0 = v0;
param.NEx = NEx;
param.m = m;
param.E = E;
param.h = h;
param.nu = nu;
param.nmode = nmode;

%% time marching
w = zeros(N/3+1,nstep); % initialize
A = zeros(nmode,nstep); % initialize
B = zeros(nmode,nstep); % initialize
w(:,1) = wini; % initial time step

for i = 2:nstep
    param.paero = paero(:,i-1);
    [~,U] = ode45(@plateDdt,[0 dt],U0,[],param);
    U0 = U(end,:)'; % retrieve solution vector U at t+dt and set it to be initial solution of the next time step
    A(:,i) = U0([1:2:end-1]); % retrieve solution A at t+dt
    B(:,i) = U0([2:2:end]); % retrieve solution B at t+dt
    w(:,i) = h .* basis(nmode,x,a) * A(:,i);
end

%% variable output
nout = max(nargout,1) - 1; % number of variable output
for i = 1:nout
    switch i
        case 1
            varargout{i} = A;
        case 2
            varargout{i} = B;
    end
end

end
