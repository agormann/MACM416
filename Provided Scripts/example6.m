clear all
%
%%%Examp
%%%Example6.m examines stability of linear multistep methods. We consider
% u'=ku,
% and u(0)=u_0. 
% The user specifies u_0, the constant k, number of time steps N
% and final time T. Final time is T=Nh.
% The output is a graph of the true solution u(t) = u_0 exp(kt), and the
% computed approximations y_n at times nh,n=0,1,2.....
% We also plot the errors at these time steps.
%% User inputs 
k=-100; % the constant in u'=ku
f=@(t,y)(k*y);

u0=1; % the initial data at t=0
T=1; % final time
dt=0.01;% number of time steps
N=T/dt; % time step size
%% Initialization
t=[0:dt:T];% creating a vector of times t_n = nh
y=zeros(N+1,1);% initializing a vector of approximations y_n
y(1)=u0; % NOTE: in Matlab, vector indices start at 1.

plot(t,exp(t),'r','LineWidth',2)

%%  A second-order accurate LMM method (method X) applied to IVP'
y(2)=exp(k*dt); % starting value Y^1
for n=2:N
    y(n+1)=-4*y(n)+5*y(n-1)+dt*( f(t(n),y(n)) + 2*f(t(n-1),y(n-1)) );
end;
y_new=y;

%% AB 2nd order

y(2)=exp(k*dt); % starting value Y^1
for n=2:N
    y(n+1)=y(n)+3/2*dt*( f(t(n),y(n)))  - (dt/2)*f(t(n-1),y(n-1)) ;
end;
y_ab=y;

%%Extrapolation from known data
y(2)=exp(k*dt); % starting value Y^1
for n=2:N
    y(n+1)=2*y(n)-y(n-1) ;
end;
y_ex=y;

%% RK4. Multistage, 4th order.
Y_RK(1)=y(1);
for n=1:N
    k1=dt*f(t(n),Y_RK(n));
    k2=dt*f(t(n)+dt/2,Y_RK(n)+k1/2);
    k3=dt*f(t(n)+dt/2,Y_RK(n)+k2/2);
    k4=dt*f(t(n)+dt,Y_RK(n)+k3);
    Y_RK(n+1)=Y_RK(n)+(k1+2*k2+2*k3+k4)/6;
end;

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
%plot(t,y_new,'ks:','LineWidth',2'),
plot(t,y_ab,'ms:','LineWidth',2'),
plot(t,y_ex,'bs:','LineWidth',2'),
plot(t,Y_RK,'gs:','LineWidth',2'),
legend('True solution u(t)','True solution at t=nh','AB 2','Extrapolated','RK4','Location','Best')
title('True and computed solutions')

% We next plot the errors
subplot(1,2,2),clc
%semilogy(t,abs(y_new'-utruepoints),'k*','LineWidth',2), hold on,
semilogy(t,abs(y_ab'-utruepoints),'m*','LineWidth',2),hold on,
semilogy(t,abs(y_ex'-utruepoints),'b*','LineWidth',2)
semilogy(t,abs(Y_RK-utruepoints),'g*','LineWidth',2)


xlabel('t'), ylabel('Error |y_n-u(nh)|')
title('Error as time progresses')
