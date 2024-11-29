% Modified p6.m script from Trefethen's 'Spectral Methods in Matlab' book
% variable coefficient wave equation with periodic boundaries
% u_t + c(x)u_x = 0 on [0,2pi]
clear;
close all;
%% Setup
nx = 256; % must be a power of 2
dx = 2*pi/nx; 
X = dx*(1:nx); 

t0 = 0;
tf = 10^4;
dt = dx/4; % size of dt vs. dx depends on method
nt = tf/dt; % this is not an integer!
nt = nx*ceil(nt/nx); % nt is an integer, and a multiple of nx (useful for plotting)
dt = tf/nt; % dt is still smaller than dx

c = 0.2 + sin(X-1).^2; % variable coefficient
k = [0:nx/2-1 0 -nx/2+1:-1]; % wave numbers

u0 = exp(-100*(X-1).^2); % initial condition
uold = exp(-100*(X-.2*dt-1).^2); % required for leapfrog method

col = [0 .5*(-1).^(1:nx-1).*cot((1:nx-1)*dx/2)]; % derivative of sinc
specD = toeplitz(col, col([1 nx:-1:2]));
F_specD = @(u) -c .* (specD*u')'; % spectral derivative (matrix)
F_fft = @(u) -c .* real(ifft(1i*k .* fft(u))); % spectral derivative (fft)
%% (a) Computing algorithm speed
t_fft = zeros(100,1);
t_specD = zeros(100,1);
for i = 1:100
    tic;
    [~,~] = main(@LF,u0,uold,dt,F_fft,nt,nx,t0);
    t_fft(i) = toc;

    tic;
    [~,~] = main(@LF,u0,uold,dt,F_specD,nt,nx,t0);
    t_specD(i) = toc;
end
t_fft_avg = mean(t_fft);
t_specD_avg = mean(t_specD);
%% Plotting solution
[U,T] = main(@RK4,u0,uold,dt,F_specD,nt,nx,t0);
%% 
close all;
figure();
view(10,70);
axis([0 2*pi T(end-nx*10) tf min(min(U)) 5]);
xlabel('x');
ylabel('t');
zlabel('u');
%colormap(1e-6*[1 1 1]);

n_to_plot = 1:nx/2:nt; % which indices to plot
hold on;
waterfall(X, T(end-nx*10:end), U(end-nx*10:end,:));
%% (b) Computing and plotting eigenvalues
ew = 1.9*nx^(-1)*eig(specD);

figure();
xlabel('Re(z)');
ylabel('Im(z)');
grid on;

hold on;
plot(real(ew), imag(ew), 'k.', 'MarkerSize', 10);
%% Time-stepping methods
function U = LF(u,uold,dt,f)
    % leapfrog
    U = uold + 2*dt .* f(u);
end

function U = RK2(u,~,dt,f)
    % runge-kutta 2nd order
    f1 = f(u);
    f2 = f(u + dt/2*f1);
    U = u + dt*f2;
end

function U = RK4(u,~,dt,f)
    % runge-kutta 4th order
    f1 = f(u);
    f2 = f(u + dt/2*f1);
    f3 = f(u + dt/2*f2);
    f4 = f(u + dt*f3);
    U = u + dt/6*(f1 + 2*(f2+f3) + f4);
end
%% Main programs
function [U,T] = main(fun,u0,uold,dt,f,nt,nx,t0)
% fun := function handle to either LF, RK2, or RK4
% u0 := initial u
% uold := time-step prior to u, only required for LF
% dt := time step
% f := function handle to either F_specD or F_fft
% nt := number of time steps to take
% nx := number of space steps to take
% t0 := initial time

U = [u0; zeros(nt,nx)];
T = [t0; zeros(nt,1)];
u = u0;
t = t0;

for i = 1:nt
    t = t+dt;
    unew = fun(u,uold,dt,f);
    uold = u;
    u = unew;
    U(i+1,:) = u;
    T(i+1) = t;
end

end
