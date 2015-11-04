function dudt = plateDdt(t,u,param)

%% parse parameters
gam = param.gam;
rho0 = param.rho;
p0 = param.p0;
a = param.a;
u0 = param.u0;
v0 = param.v0;
NEx = param.NEx;
m = param.m;
E = param.E;
h = param.h;
nu = param.nu;
nmode = param.nmode;
paero = param.paero;

%% nondimensionalize
D = E*h^3 / (12*(1-nu^2));
Rx = NEx*a^2 / D;

%% evaluate dudt
A = u([1:2:end-1]); % [A1,A2,...]'
B = u([2:2:end]); % [A1',A2',...]'
modeId = [1:nmode]'; % mode index; [1,2,...,nmode]'
Ar_sum = sum(A.^2 .* modeId.^2 * pi^2 / 2); % precompute

dAdt = B;

intPaero = integratePaero(paero,a,nmode,modeId);
dBdt = ( -(modeId*pi).^4 - 6*(1-nu^2)*Ar_sum*(modeId*pi).^2 ...
    - Rx*(modeId*pi).^2 ) .* A - 2*a^3/(D*h)*intPaero';

dudt = zeros(2*nmode,1);
dudt([1:2:end-1]) = dAdt;
dudt([2:2:end]) = dBdt;

end