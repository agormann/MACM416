% example16: spatial differences. We compare 4 different methods for the
% second derivative
clear all, clc,clf
%% Problem specification
% x in [-M,M], x_0=-M, x_{N}=M
% Due to periodic boundary conditions we don't need x_0
N=200; 
h=2*pi/(N); xgrid=h*(1:N)';x=xgrid; 
M=pi
%test  data for differentiation.

phi=@(x) (-cos(x)); true=@(x)(cos(x));
%phi=@(x)(cos(2*x)); true=@(x)(-4*cos(2*x));

figure(1),hold on

%plot(xgrid,phi(xgrid),'ko','LineWidth',2),hold on

%%
% Spatial differencing on periodic grid
% For single derivative: u_x(x_i,t) replaced by (u(x_i+1,t)-u(x_i-1,t) )/2h
firstcol=[0 -1 zeros(1,N-2)]'; firstrow=[0 1 zeros(1,N-2)];
SD=toeplitz(firstcol,firstrow);
SD(1,N)=-1;SD(N,1)=1; % modifications to make differencing periodic
SD=SD/(2*h);
SD2=SD*SD;% second derivative using repeated first derivatives

err=SD2*phi(x)-true(x);
plot(x,abs(err),'b','LineWidth',2)


 %% for second order derivative using 2-point centered scheme
firstcol=[-2 1 zeros(1,N-2)]'; firstrow=[-2 1 zeros(1,N-2)];
TD=toeplitz(firstcol,firstrow);
TD(1,N)=1;TD(N,1)=1; % modifications to make differencing period
TD=TD/h/h;
err=TD*phi(x)-true(x);
plot(x,abs(err),'g','LineWidth',2)

    %% Construct spectral differentiation matrix from explicit derivative matrix 
  tstart=tic;
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
    SpecD = toeplitz(column,column([1 N:-1:2]));
    SpecD=SpecD*SpecD;
    W=SpecD*phi(x);
    
    telapsed_explicit=toc(tstart)
    
    err=W-true(x);
    
    plot(x,abs(err),'r')
 
    
    
    %% Construct spectral differentiation matrix using FFT
 
tstart=tic;
U=fft(phi(x));    
W1=  (ifft(1i*[0:N/2-1 0 1-N/2:-1]'.*U)); % diff wrt x
W2 = (ifft(-[0:N/2 1-N/2:-1]'.^2.*U)); %second diff wrt x

telapsed_fft=toc(tstart)

err=W2-true(x);
 plot(x,abs(err),'co-')

   %%
   figure(1), xlabel('x'), ylabel('Absolute error in 2nd derivative'), legend('Using SD2', 'Using TD', 'Using spectral 1', 'Using spectral FFT')
    title('Checking application of differentiation matrix')
    
