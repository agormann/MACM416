clear all
%
%%%Examp
%%%Example1.m is an implementation of the Forward Euler method, along with
% inline function specifications, to solve (IVP')
% u'=ku,
% and u(0)=u_0. 
% The user specifies u_0, the constant k, number of time steps N
% and final time T. Final time is T=Nh.
% The output is a graph of the true solution u(t) = u_0 exp(kt), and the
% computed approximations y_n at times nh,n=0,1,2.....
% We also plot the errors at these time steps.
%% User inputs 
k=-100; % the constant in u'=ku
u0=1; % the initial data at t=0
T=2; % final time
N=100;% number of time steps
h=T/N; % time step size
%% The Forward Euler method applied to IVP'
t=[0:h:T];% creating a vector of times t_n = nh
y=zeros(N+1,1);% initializing a vector of approximations y_n
y(1)=u0; % NOTE: in Matlab, vector indices start at 1.
t(1)=0 ;
for n=1:N
    y(n+1)=y(n)+k*h*y(n);% the Euler update
end

%% Plots
% we first plot the true solution, and explicitly mark the solution at the
% times t_n:
s=linspace(0,T,1e6); utrue=u0*exp(k*s);utruepoints=u0*exp(k*t);
figure(1), clf,

subplot(1,2,1), 
plot(s,utrue,'r','LineWidth',2'),hold on
xlabel('t'),ylabel('y')
plot(t,utruepoints,'bo','LineWidth',2'),
%%

% Next we overlay the approximated solution, and join the points by
% straight line segments
plot(t,y,'ms:','LineWidth',2'),
legend('True solution u(t)','True solution at t=nh','Computed solution y_n','Location','Best')
title('True and computed solutions')

% We next plot the errors
subplot(1,2,2),clc
semilogy(t,abs(y'-utruepoints),'k*','LineWidth',2') 
xlabel('t'), ylabel('Error |y_n-u(nh)|')
title('Error as time progresses')