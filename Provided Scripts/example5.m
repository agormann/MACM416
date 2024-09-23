%%example 5: drawing regions of stability for initial value solvers


%% Runge-Kutta methods: if Y^n+1 = R(z)Y^n, find z so that |R(z)|<=1.
N=100;
x = linspace(-4,4,N);y = linspace(-4,4,N);

[X,Y] = meshgrid(x,y);
Z = X + i*Y; % generating a grid in the complex plane

% 2nd order Runge-Kutta method
RK2 = @(z) (1 + z + 1/2*z.^2);

% 3rd order Runge-Kutta method
RK3 = @(z) (1 + z + 1/2*z.^2 + 1/6*z.^3);

% 4th order Runge-Kutta method
RK4 = @(z) (1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4);

% Evaluation of R(z) 
Rval2 = abs(RK2(Z)); Rval3 = abs(RK3(Z)); Rval4 = abs(RK4(Z));

figure(1),surfc(X,Y,Rval3), xlabel('real'), ylabel('imag'),colorbar, shading flat


%%
figure(2), clc, 
plot([-5,5 ],[0 0],'k-','LineWidth',2), hold on,
plot([0 0],[-5 5],'k-','LineWidth',2)
contour(x,y,Rval2,[1 1],'r-','LineWidth',2); 
contour(x,y,Rval3,[1 1],'b-','LineWidth',2);
contour(x,y,Rval4,[1 1],'g-','LineWidth',2);
legend('RK2','RK3','RK4')
title('Region of absolute stability for RK methods','FontSize',15)

%% Multistep methods. This code is taken from 
% p25.m, found in 'Spectral Methods in Matlab' by L. N Trefethen
% p25.m - stability regions for ODE formulas

figure(3)
% Adams-Bashforth:
  clf, subplot('position',[.1 .56 .38 .38])
%  plot([-8 8],[0 0]), hold on, plot([0 0],[-8 8])
  z = exp(1i*pi*(0:200)/100); r = z-1;
  s = 1; plot(r./s),             hold on                     % order 1
  s = (3-1./z)/2; plot(r./s)                         % order 2
  s = (23-16./z+5./z.^2)/12; plot(r./s)              % order 3
  axis([-2.5 .5 -1.5 1.5]), axis square, grid on
  title Adams-Bashforth
legend('Order 1','Order 2', 'Order 3');
%%
% Adams-Moulton:
  subplot('position',[.5 .56 .38 .38])
  % plot([-8 8],[0 0]), hold on, plot([0 0],[-8 8])
  s = (5*z+8-1./z)/12; plot(r./s), hold on                    % order 3
  s = (9*z+19-5./z+1./z.^2)/24; plot(r./s)           % order 4 
  s = (251*z+646-264./z+106./z.^2-19./z.^3)/720; plot(r./s)    % 5
  d = 1-1./z;
  s = 1-d/2-d.^2/12-d.^3/24-19*d.^4/720-3*d.^5/160; plot(d./s) % 6
  axis([-7 1 -4 4]), axis square, grid on, title Adams-Moulton
legend('Order 3', 'Order 4', 'Order 5', 'Order 6')
