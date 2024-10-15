% example12: a first look at the method of lines
% For each of the following examples, we assume periodic boundary
% conditions. For comparison, we'll use the same  initial condition
% example 1: 1st order wave equation u_t= cu_x
% example 2: Heat equation, u_t = cu_xx
% example 3: Schrodinger equation, u_t = icu_xx
clear; clc; clf;
%% Problem specification
% x in [-M,M], x_0=-M, x_{N+1}=M
% Due to periodic boundary conditions we don't need x_0
M=10; N=60; h=2*M/(N+1); xgrid=[-M+h:h:M];

% c=speed
c=0.5;

%Initial data.
phi=@(x)sin(pi*x/M)+0.25*cos(5*pi*x/M);
% Initial data discontinuous:For fun, we'll try a 'hat'
 % piece1=@(x) abs(x)<=1;
  %piece2=@(x) abs(x)>1;
%  
 % phi= @(x) piece1(x) + piece2(x).*0;


figure(1)

plot(xgrid,phi(xgrid),'LineWidth',2),hold on

%%
% Spatial differencing on periodic grid
% For single derivative: u_x(x_i,t) replaced by (u(x_i+1,t)-u(x_i-1,t) )/2h
firstcol=[0 -1 zeros(1,N-1)]'; firstrow=[0 1 zeros(1,N-1)];
SD=toeplitz(firstcol,firstrow);
SD(1,N+1)=-1;SD(N+1,1)=1; % modifications to make differencing periodic
SD=SD/(2*h);

 % for second order derivative
firstcol=[-2 1 zeros(1,N-1)]'; firstrow=[-2 1 zeros(1,N-1)];
TD=toeplitz(firstcol,firstrow);
TD(1,N+1)=1;TD(N+1,1)=1; % modifications to make differencing periodic
TD=TD/h/h;
%% Wave equation u_t = cu_x by the method of lines. 
InitialVector=phi(xgrid)';
wave=@(t,Y)(c*SD*Y);
solwave=ode23(@(t,y)wave(t,y),[0,20],InitialVector);
figure(2), hold on,
for j=0:0
plot(xgrid,deval(solwave,j),'LineWidth',2),
end
legend('t=0','t=1','t=2','t=3','t=4','t=5')
title('Linear Advection'), xlabel('x')
%% Heat equation

heat=@(t,Y)(c*TD*Y);
solheat=ode45(@(t,y)heat(t,y),[0,20],InitialVector);
figure(3), hold on,
for j=0:5
plot(xgrid,deval(solheat,j),'LineWidth',2),
end
legend('t=0','t=1','t=2','t=3','t=4','t=5')
title('Linear heat equation'), xlabel('x')

%% Schrodingers equation
schrodinger=@(t,Y)(1i*c*TD*Y);
solschrodinger=ode45(@(t,y)schrodinger(t,y),[0,20],InitialVector);
figure(4), subplot(2,1,1),hold on,subplot(2,1,2),hold on,
for j=0:5
subplot(2,1,1)
plot(xgrid,real(deval(solschrodinger,j)),'LineWidth',2),
subplot(2,1,2)
plot(xgrid,imag(deval(solschrodinger,j)),'LineWidth',2),
end
subplot(2,1,1)
legend('t=0','t=1','t=2','t=3','t=4','t=5')
title('Real part: Linear Schrodinger equation'), xlabel('x')
subplot(2,1,2)
legend('t=0','t=1','t=2','t=3','t=4','t=5')
title('Imaginary part: Linear Schrodinger equation'), xlabel('x')
