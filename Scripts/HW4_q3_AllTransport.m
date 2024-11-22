% Solves the linear advection equation over the interval [0,1] using 
% one of three different schemes: 'Lax-Friedrichs', 'upwind', or 'Lax-Wendroff'. 
% Periodic boundary conditions are assumed.
%% Parameters
%   a - Advection coefficient (velocity)
%   nu - Local Courant number. Determines Courant number and temporal grid size.
%   theta - The power which nu is raised to. Determines whether the scheme
%     used is 'Lax-Friedrichs' (-1), 'upwind' (0), or 'Lax-Wendroff' (1).
%   N - Number of grid points. Determines spatial grid size dx.
%   NumIts - Number of iterations of the scheme.
%   key - Chooses the initial condition. It is either the indicator function of the 
%     interval [1/4,3/4] (1), the hat function of the interval [1/4,3/4] (2), 
%     or sin(2*pi*x) (otherwise).

a = 1;
nu = 1.01;
theta = 1;
N = 100;
NumIts = 25;
key = 2;
%% Setup
h = 1/N; % spatial grid size
x = h/2 : h : 1-h/2; % discretization of the interval [0,1]
lambda = nu/abs(a); % CFL number
dt = lambda*h; % temporal grid size

% Initial condition
if key == 1
    u = 0 + ((1/4-abs(x-1/2))>0);
elseif key == 2
    u = max(0, 1-4*abs(x-1/2));
else
    u = sin(2*pi*x);
end
plot(x, u, '.');
%% Numerical Solution
time = 0;
for it = 1:NumIts
    time = time + dt;
    % computing numerical flux function
    Phi = 0.5*a*([u, u(1)] + [u(N), u]) ...
    - 0.5*nu^theta*abs(a)*([u, u(1)] - [u(N), u]);
    % updating u with forward time-step
    u = u - lambda * (Phi(2:N+1) - Phi(1:N));
    
    xmodulo = mod(a*time,1);
    % computing exact solution for error comparison
    if key == 1
        uex = 0 + ((1/4-abs(x-xmodulo-1/2))>0) ...
        + ((1/4-abs(x-xmodulo+1/2))>0);
    elseif key == 2
        uex = max(0, 1-4*abs(x-xmodulo-1/2)) ...
        + max(0, 1-4*abs(x-xmodulo+1/2));
    else
        uex = sin(2*pi*(x-xmodulo));
    end

    if mod(it, 50)
        subplot(1,2,1), plot(x, u, '.', x, uex, '-');
        title('Discrete solution (uh)');
        subplot(1,2,2), plot(x, u-uex, '.');
        title('Error (uh-u_ex)');
        
        drawnow();
    end
end