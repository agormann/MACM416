function [xx,y,h,D] = poisson1D(x,n,f,BC)
%POISSONODE Solves the ODE y''(x) = f(x) with boundary conditions y(a) = alpha, y'(b) = beta.
%   PARAMETERS:
%   x  : [a b]
%   n  : dimension of D
%   f  : function handle
%   BC : [alpha beta]
% 
%   SOLUTION:
%   xx : interpolation of [a,b]
%   y  : solution to ODE
%
%   EXPLANATION:
%   Rewrite problem as D(Y)+B = F, where
%   D : 3-point centered differentiation matrix
%   y : approximation to y(x_n)
%   B : boundary conditions
%   F : value of f(x_n)

% setup
h = (x(2)-x(1)) / (n-1);
xx = (x(1):h:x(2))';
alpha = BC(1);
beta = BC(2);

% spatial differencing
D = spdiags([1 -2 1], [-1 0 1], n-1, n-1);

% boundary conditions
D(end,:) = [zeros(1,n-3) 2 -2];
B = zeros(n-1,1);
B(1) = alpha/h^2;
B(end) = 2*beta/h;

% solving
if isa(f, 'function_handle')
    F = f(xx(2:end)); % x=a is pinned to y=alpha
    D = D/h^2;
    y = D\(F-B); % dimension (n-1)x1
    y = [0; y]; % needs to be nx1
else
    D = D/h^2;
    y = null;
end

end
