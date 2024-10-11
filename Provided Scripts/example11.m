%% example 11: a first look at BVP

%% we'll be solving u''=f(x), u(a)=alpha, u(b)=beta
clear all, clc
% problem parameters

f=@(x)(sin(x));
alpha=0;beta=pi;
utrue=@(x)(-sin(x)+x); % true solution for this problem

%% Set up the grid a=x_0<x_1<...x_N<x_N+1=b
a=0; b=pi;
figure(1), hold on, xlabel('x'), title('True solution in blue')
figure(2)

%Numpoints=5;  
Numpoints=5*2.^[1:4];

for i=1:length(Numpoints)
   N=Numpoints(i); 
   h=(b-a)/(N+1);% grid size
   xgrid=[a+h:h:b-h]';% grid points x1,x2,...x_N
% setting up matrix and right hand side
   c=[-2 1 zeros(1,N-2)]; A=(1/h^2)*toeplitz(c); % matrix for centered difference of d^2/dx^2
   Fh=zeros(1,N);
   Fh=f(xgrid); Fh(1)=Fh(1)-alpha/h^2; Fh(end)=Fh(end)-beta/h^2;

%% Solving the problem
Ugrid=A\Fh; 
xfull=[a; xgrid; b];% adding on end points to our grid
Ufullgrid=[alpha; Ugrid; beta]; utruefullgrid=utrue(xfull);
figure(1)
%plot(xfull, Ufullgrid,'o','LineWidth',2,'MarkerFaceColor','r'); 
plot(eig(A),'r*')
xfine=linspace(a,b,1000);
%plot(xfine, utrue(xfine),'b','LineWidth',2); hold on;

Errorvector=Ufullgrid-utruefullgrid;

err1(i)=norm(Errorvector,1); err2(i)=norm(Errorvector,2);
griderr1(i)=h*err1(i); griderr2(i)=sqrt(h)*err1(i); % grid norms 


end;
%%

%loglog(Numpoints,err1,'r*',Numpoints,griderr1,'rs',Numpoints,err2,'b*',Numpoints,griderr2,'bs')
xlabel('number of points'), ylabel('error'), legend('L1 error', 'Grid-L1 error','L2 error','Grid-L2 error')