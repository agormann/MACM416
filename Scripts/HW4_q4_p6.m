% Modified p6.m Script
% variable coefficient wave equation from Trefethen's 'Spectral Methods in Matlab' book
% u_t + c(x)u_x = 0 with periodic boundaries
clear;
close all;
%% 
N = 128;
h = 2*pi/N; 
x = h*(1:N); 
t = 0; 
dt = h/4;
c = 0.2 + sin(x-1).^2; % variable coefficient
v = exp(-100*(x-1).^2); % initial condition
vold = exp(-100*(x-.2*dt-1).^2);
k = [0:N/2-1 0 -N/2+1:-1];
%% 
tmax = 20;
tplot = .50;
plotgap = round(tplot/dt); 
dt = tplot/plotgap;
nplots = round(tmax/tplot);
data = [v; zeros(nplots,N)]; 
tdata = t;
U = [v; zeros(nplots*plotgap-1,N)];
T = [t; zeros(nplots*plotgap-1,1)];

for i = 1:nplots
    for n = 1:plotgap
        t = t+dt;
        w_hat = 1i*k .* fft(v);
        w = real(ifft(w_hat));
        vnew = vold - 2*dt*c.*w; 
        vold = v; 
        v = vnew;
        U(plotgap*(i-1)+n,:) = v;
        T(plotgap*(i-1)+n) = t;
    end
    data(i+1,:) = v; 
    tdata = [tdata; t];
end
%%
close all;
figure();
view(10,70);
view(3)
axis([0 2*pi 0 tmax min(min(data)) 5]);
ylabel('t');
zlabel('u');
%colormap(1e-6*[1 1 1]);

hold on;
%waterfall(x,tdata,data);
waterfall(x,T,U);