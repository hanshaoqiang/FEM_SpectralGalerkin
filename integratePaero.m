function intPaero = integratePaero(paero,a,nmode,modeId)
% This function computes the integral of \Delta p_aero * phi_modeId over
% [0,a]. Note that the integral is evaluated approximately.
%% input
% paero: \Delta p_aero
% a: plate length
% nmode: number of modes
% modeId: denote the index of basis function in the integrand; can be a
% scalar or column vector

%% output
% intPaero: size = [length(modeId),1] row vector. Each entry is the
% integral corresponding to the basis function phi_modeId.

%% execution
N = length(paero) - 1; % number of horizontal (x direction) elements across the flow field
x = linspace(0,a,N/3+1)'; % coordinates of nodes on the plate
if abs( max(diff(x)) - min(diff(x)) ) < 10*eps % sanity check
    dx = mean(diff(x));
else
    error('dx is not uniform.')
end

if 0
    paero_plate = paero([N/3+1:2*N/3+1]); % pick out \Delta p_aero over the plate
    paero_plate_avg = 0.5.*(paero_plate(1:end-1) + paero_plate(2:end)); % taking average paero each cell
    x_avg = 0.5 .* (x(1:end-1) + x(2:end)); % coordinates of cell center
    Phi_avg = basis(nmode,x_avg,a); % evaluate basis functions at x_avg
    
    % approximate the integral by integrating piecewise constant approximation
    intPaero = sum(repmat(paero_plate_avg,1,length(modeId)) .* Phi_avg(:,modeId) .* dx , 1);
else
    % integrate analytically: only applicable for the case in the project
    % description
    intPaero = -paero(1) * a ./ (modeId'*pi) .*( (-1).^(modeId') - 1 );
end

end