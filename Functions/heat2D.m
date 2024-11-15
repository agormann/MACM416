function [t,U,X,Y] = heat2D(x,N,T,options)
%heat2D Solves the heat equation in 2D
%   PARAMETERS:
%       x : [x0 xend]
%           assuming square region! 
%       N : number of points
%       T : [t0 tend]
%       options : ['z'||'p', 'z'||'p' 'z'||'p'] 
%                 chooses source term, initial conditions, and boundary
%                 conditions. 'z' := zeros and 'p' := peaks.
%   SOLUTION:
%       t : column vector of times
%       U : vertically concatenated solutions in the form of row vectors
%       X : mesh for plotting
%       Y : mesh for plotting
%   EXPLANATION
%

% setup
n = N-2;
h = (x(2)-x(1)) / (N-1);
xx = (x(1):h:x(2))';
[X,Y] = meshgrid(xx,xx);

% spatial differencing matrix
D = -delsq(numgrid('S',N));

% solving for source
if options(1) == "z"
    F = sparse(n^2,1);
elseif options(1) == "p"
    F = peaks(X(2:end-1, 2:end-1), Y(2:end-1, 2:end-1));
    F = reshape(F,n^2,1);
else
    error("Parameter source must be either zero or peaks.")
end

% solving for initial
if options(2) == "z"
    U0 = zeros(n^2,1);
elseif options(2) == "p"
    U0 = peaks(X(2:end-1, 2:end-1), Y(2:end-1, 2:end-1));
else
    error("Parameter IC must be either zero or peaks.")
end

% solving for boundary
if options(3) == "z"
    B = sparse(n^2,1);
    Z = sparse(N,N);
elseif options(3) == "p"
    B = sparse(n^2,N^2);
    Z = peaks(X,Y);

    % painstakingly found indices :)
    i = 1:n;
    idFirst = i*(n^2+1);
    idLast = ((N-1)*N+i)*n^2 + (n-1)*n + i;
    idLeading = i*N*n^2 + (i-1)*n + 1;
    idTrailing = ((i+1)*N-1)*n^2 + i*n;
    
    id = [idFirst idLast idLeading idTrailing];
    B(id) = 1; % do spy(B) for structure

    % compute *full* X,Y then stack it vertically
    B = B*reshape(peaks(X,Y),N^2,1);
else
    error("Parameter BC must be either zero or peaks.")
end

% solution
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t,y] = ode45(@(t,y) ode(t,y,D,F,B), T, U0, options);

% including boundary in solution
U = zeros(length(t),N^2);
for i = 1:length(t)
    u = y(i,:);
    u = reshape(u,n,n);
    Z = reshape(Z,N,N);
    Z(2:end-1,2:end-1) = u;
    Z = reshape(Z,1,N^2);
    U(i,:) = Z;
end

end

% required for ode45
function dydt = ode(t,y,D,F,B)
    dydt = D*y-F+B;
end