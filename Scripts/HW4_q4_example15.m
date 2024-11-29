%% example 15: wave equation in 2 dimension
% we solve u_t = c_1 u_x + c_2 u_y with periodic boundary conditions
% the bulk of this code is about setting up difference operators for the 
% periodic boundary condition case on [-A,A] X [-B,B]


% problem inputs
tmax=5;
c1= 3; c2=0; % wave speed
%phi=@(x,y) exp(-20*(x).^2-20*y.^2)+0.5*exp(-20*(x-pi/2).^2 -20*(y-pi/4).^2 ); % initial data 
phi=@(x,y) exp(-20*(x.^2+y.^2)); % initial data

A=pi;Nx=15; % for x direction
B=pi; Ny=15; % for y direction
dt=0.01; 

nmax= round(tmax/dt); % total number of time steps
nplt=floor((tmax/8)/dt); % we'll plot 20 time snapshots of the solutions

%% % setting up the spatial grid and differentiation in x direction on [-A:A], x_0=-A, x_{Nx+1}= A

dx=2*A/(Nx+1); x=[-A+dx:dx:A]';

% Dx= centered differencing for first derivative
firstcol=[0 -1 zeros(1,Nx-1)]'; firstrow=[0 1 zeros(1,Nx-1)];
Dx=toeplitz(firstcol,firstrow);
Dx(1,Nx+1)=-1;Dx(Nx+1,1)=1; % modifications to make differencing periodic
Dx=Dx/(2*dx);
Ix=eye(Nx+1);



%% % setting up the spatial grid in y direction on [-B:B], y_0=-B, y_{Ny+1}= B
 
dy=2*B/(Ny+1); y=[-B+dy:dy:B]';
% Dy= centered differencing for first derivative
firstcol=[0 -1 zeros(1,Ny-1)]'; firstrow=[0 1 zeros(1,Ny-1)];
Dy=toeplitz(firstcol,firstrow);
Dy(1,Ny+1)=-1;Dy(Ny+1,1)=1; % modifications to make differencing periodic
Dy=Dy/(2*dy);
Iy=eye(Ny+1);
%%
%%% Combining onto tensorial grid in 2 dimensions

delta_x = kron(Dx,Iy); % this is the x derivative in a 2-dimensional setting
delta_y = kron(Ix,Dy); % this is the y derivative in a 2-dimensional setting.

 [XX,YY] = meshgrid(-A:A/100:A,-B:B/100:B); % finer grid for plotting
%% testing the derivatives
[X,Y]=meshgrid(x,y); % tensorial grid
Utest=sin(X).*cos(2*Y); % Utest is a matrix the same size as X or Y

Utest_x = cos(X).*cos(2*Y); Utest_y=-2*sin(X).*sin(2*Y); % true derivatives

% test of x derivative
U_vec = Utest(:); % reshaping the grid function into a vector
U_x_num = delta_x*U_vec; % Apply the derivative matrix in x direction
U_x_matrix = reshape(U_x_num,Ny+1,Nx+1); % reshaping vector back into grid function

%%
figure(1)
subplot(1,2,1)
contourf(X,Y,Utest_x)
colorbar
title 'Numeric U_{x}'
subplot(1,2,2)
contourf(X,Y,abs(Utest_x-U_x_matrix))
colorbar
title 'error with respect to Analytic U_{x}'

%% 
% test of y derivative

U_vec = Utest(:);
U_y_num = delta_y*U_vec;

U_y_matrix = reshape(U_y_num,Ny+1,Nx+1);

%%
figure(2)
subplot(1,2,1)
contourf(X,Y,U_y_matrix)
colorbar
title 'Numeric U_{y}'
subplot(1,2,2)
contourf(X,Y,abs(Utest_y-U_y_matrix))
colorbar
title 'error with respect to Analytic U_{y}'

%%%
%%% solving the PDE u_x = c_1 u_x + c_2 u_y, u(x,y,0)= phi(x,y)
%% collecting data for plotting
tdata=0;

initdata=phi(X,Y); initdata=initdata(:);

truedata=initdata; % true solution

LWn=initdata; % Forward euler
figure(3)
for n=1:nmax
    t=n*dt;
    LWn=LWn+dt*[c1*delta_x + c2*delta_y]*LWn; % Forward Euler step
    if mod(n,nplt)==0
          
         plotsnap=reshape(LWn,Ny+1,Nx+1);
        
      vvv = interp2(X,Y,plotsnap,XX,YY,'cubic');
      mesh(XX,YY,vvv),
      view([10,70])
        % mesh(X,Y,plotsnap);
        xlabel x, ylabel y, 
        drawnow
         tdata=[tdata t];
        pause
    end
  
end


