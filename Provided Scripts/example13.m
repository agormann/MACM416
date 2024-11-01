% example13: spatial difference
%clear all, clc,clf
%% Problem specification
% x in [-M,M], x_0=-M, x_{N}=M
% Due to periodic boundary conditions we don't need x_0
N=12; 
h=2*pi/(N); xgrid=[-pi+h:h:pi];
M=pi
%test  data for differentiation.

phi=@(x) (-cos(x)); true=@(x)(cos(x));
% phi=@(x)(cos(2*x)); true=@(x)(-4*cos(2*x));

figure(1),hold on

%plot(xgrid,phi(xgrid),'ko','LineWidth',2),hold on

%%
% Spatial differencing on periodic grid
% For single derivative: u_x(x_i,t) replaced by (u(x_i+1,t)-u(x_i-1,t) )/2h
h = 2*pi/N; x = -pi + (1:N)'*h; % (A) WHAT IS HAPPENING HERE?
firstcol=[0 -1 zeros(1,N-2)]'; firstrow=[0 1 zeros(1,N-2)];

SD=toeplitz(firstcol,firstrow);
SD(1,N)=-1;SD(N,1)=1; % modifications to make differencing periodic
SD=SD/(2*h);
SD2=SD*SD;% second derivative using repeated first derivatives

err=SD2*phi(x)-true(x);
plot(x,abs(err),'b','LineWidth',2)


 %% for second order derivative
firstcol=[-2 1 zeros(1,N-2)]'; firstrow=[-2 1 zeros(1,N-2)];
TD=toeplitz(firstcol,firstrow);
TD(1,N)=1;TD(N,1)=1; % modifications to make differencing periodic
TD=TD/h/h;
err=TD*phi(x)-true(x);
plot(x,abs(err),'g','LineWidth',2)

    %% Construct spectral differentiation matrix:
    h = 2*pi/N;
    x = -pi + (1:N)'*h;
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
    SpecD = toeplitz(column,column([1 N:-1:2]));
    SpecD=SpecD*SpecD;
    err=SpecD*phi(x)-true(x);
    plot(x,abs(err),'r')
 
    figure(1), xlabel('x'), ylabel('Absolute error in 2nd derivative'), legend('Using SD2', 'Using TD', 'Using spectral')
    title('Checking application of differentiation matrix')
    %% Test on Poisson problem u''= f with periodic BC
    
    f=@(x)sin(x); % true solution is -sin(x) + C
    % f=@exp(-(x.^2)/100);
    figure(2), hold on
   
    usd2=SD2\f(x);  plot(x,usd2,'b','LineWidth',2);% using SD2  
    ut=TD\f(x);  plot(x,ut,'g','LineWidth',2);% using TD
    uspec=SpecD\f(x);  plot(x,uspec,'r','LineWidth',2);% using Spectral
   
    figure(2), xlabel('x'), ylabel('Solution of PDE'), legend('Using SD2', 'Using TD', 'Using spectral')
   title('Solving Poisson problem')