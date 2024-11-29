% Parameters
N = 128;             % Number of spatial grid points
h = 2*pi/N;            % Grid spacing
x = h*(1:N); 
dt = h/4;
t_final = 5;         % Final time

% Wavenumbers for FFT
k = [0:N/2-1 0 -N/2+1:-1];

% Initial condition
u0 = exp(-100*(x-1).^2);

% Time-stepping setup
u = u0;
t = 0;

% Preallocate storage for visualization
plot_every = 50; % Plot every 50 time steps
time_steps = floor(t_final/dt);
u_storage = zeros(length(x), ceil(time_steps/plot_every) + 1);
u_storage(:,1) = u0;

% Define RHS of the wave equation in spectral space
function du_dt = rhs(u_hat, k, x)
    u = ifft(u_hat);               % Transform back to physical space
    u_x_hat = 1i * k .* u_hat;     % Spectral derivative
    u_x = ifft(u_x_hat);           % Transform derivative back
    du_dt = -fft(x'.*u_x);         % x*u_x term in spectral space
end

% Time stepping
step = 1;
for n = 1:time_steps
    u_hat = fft(u);  % Compute FFT of solution
    
    % Runge-Kutta 4th order scheme
    k1 = dt * rhs(u_hat, k, x);
    k2 = dt * rhs(u_hat + 0.5*k1, k, x);
    k3 = dt * rhs(u_hat + 0.5*k2, k, x);
    k4 = dt * rhs(u_hat + k3, k, x);
    u_hat_new = u_hat + (k1 + 2*k2 + 2*k3 + k4)/6;

    % Transform back to physical space
    u = ifft(u_hat_new, 'symmetric');
    t = t + dt;
    
    % Store solution for visualization
    if mod(n, plot_every) == 0
        step = step + 1;
        u_storage(:, step) = u;
    end
end

% Visualization
[X, T] = meshgrid(x, linspace(0, t_final, size(u_storage, 2)));
figure;
surf(X, T, u_storage', 'EdgeColor', 'none');
view(2);
xlabel('x');
ylabel('t');
title('Wave Equation Solution');
colorbar;
