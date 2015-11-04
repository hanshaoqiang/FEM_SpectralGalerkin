function phi = basis(nmode,x,a)
% This function computes nmode Fourier basis function values at coordinates
% x, where x can be a scalar or column vector

% x must be a scalar or column vector
if size(x,2) > 1
    error('x is not a column vector.')
end

% initialize
phi = zeros(size(x,1),nmode);

% evaluate basis functions
for i = 1:nmode
    phi(:,i) = sin(i*pi*x/a);
end