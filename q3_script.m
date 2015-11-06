close all; clear all

endT = 3;
nstep = 31;

nmode = 5;
N = 30;

%% test input
rho0 = 1; p0 = 1/1.4; a = 2;
u0 = 2; 
NEx = -0.526;
m = 100;
E = 72.8E6;
h = 0.002;
nu = 0.3;

if (N/3+1) < 2*nmode
    error('too many modes for the mesh resolution. aliasing might happen')
end
wini = zeros(N/3+1,1);
wtini = zeros(N/3+1,1);

t = linspace(0,endT,nstep)';
paero = repmat( (t>=0 & t<=10)', N+1,1) ...
    - repmat( (t>10 & t<=30)', N+1,1);

dx = 3*a/N;
CFL = u0*endT / ((nstep-1)*dx)
%%
[w,A,B] = q4(rho0, p0, a, u0, NEx, m, E, h, nu,...
    endT, nstep, paero, nmode, wini, wtini);

%% post processing
% plot w at t 0, 10, 20, 30
x = linspace(0,a,N/3+1)';
phi_x = basis(nmode,x,a);
xx = linspace(0,a)';
phi_xx = basis(nmode,xx,a);

dNt = (nstep - 1)/3;

% % compare ways of presenting w
% for i = 1:3
%     w_discrete = w(:,dNt*i+1);
%     
%     A_approx = phi_x' * w_discrete * dx*2/(a*h);
%     
%     w_series = h * phi_xx * A(:,dNt*i+1);
%     w_approx = h * phi_xx * A_approx;
%     
%     figure(i)
%     plot(x,w_discrete,'k+-', xx,w_approx,'rx', xx,w_series,'bo')
% end

% % t = 0
% figure()
% plot(x,wini,'r+-')

ctype = {'-k', '-b', '-r', '-y', 'x-c', 's-g', 'p-y', '>-k'};

% w for t = 0, 10, 20, 30
figure(1)
set(gca,'FontSize',12)
for i = 1:4
    if i == 1
        w_series = interp1(x,wini,xx,'linear');
    else
        w_series = h * phi_xx * A(:,dNt*(i-1)+1); 
    end    
    plot(xx,w_series,ctype{i})
    hold on
end
legend({'$t = 0$','$t = 10$','$t = 20$','$t = 30$'},...
    'Location','best','Interpreter','latex')
xlabel('$x$','interpreter','latex')
ylabel('$w$','interpreter','latex')

fname = sprintf('w.eps');
print('-depsc2',fname);
unix(sprintf('epstopdf %s', fname));
delete(fname); % delete eps files

% solution w & wt @ x=1.5 for t=[0,20]
x0 = 1.5;
phi_x0 = basis(nmode,x0,a);
w_x0 = h * phi_x0 * A(:,[1:2*dNt+1]);

D = E*h^3 / (12*(1-nu^2));
dtaudt = sqrt( D / (m*a^4) ); % dtau/dt
wt_x0 = h * phi_x0 * dtaudt * B(:,[1:2*dNt+1]);

dt = endT / (nstep-1);

% figure()
% plot([0:dt:2*dNt*dt],w_x0,'r+-')
% 
% figure()
% plot([0:dt:2*dNt*dt],wt_x0,'r+-')

figure(2)
set(gca,'FontSize',12)
plot(w_x0,wt_x0,'k+-')
xlabel('$w\vert_{x=1.5}$','interpreter','latex')
ylabel('$\frac{\partial w}{\partial t}\vert_{x=1.5}$','interpreter','latex')

fname = sprintf('dwdt_w.eps');
print('-depsc2',fname);
unix(sprintf('epstopdf %s', fname));
delete(fname); % delete eps files