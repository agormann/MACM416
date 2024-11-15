%%% For comparing three different schemes for the transport equation 
%u_t +a u_x =0.

%%
a=-1;

theta = 1;
N = 1000;
h = 1 / N;
x = [h/2 : h: 1-h/2];
key =2; % 1 for an initial condition which is a step function
NumIts=200;
nu=1.01;
time=0;
 
 
 %%
 
if (key==1)
u = 0 + ((1/4-abs(x-1/2))>0); % IC is a step function
elseif (key==2)
 u = max(0, 1-4*abs(x-1/2));
 else
 u = sin(2*pi*x);
end;
plot(x, u, '.');


lambda = nu / abs(a);
dt = h * lambda;

for it=1:NumIts
time = time + dt;
Phi = 0.5*a*([u, u(1)] + [u(N), u]) ...
- 0.5*nu^theta*abs(a)*([u, u(1)] - [u(N), u]);

u = u - lambda * (Phi(2:N+1) - Phi(1:N));

 xmodulo = mod(a*time,1);
 % exact solutions
if (key==1) 
 uex = 0 + ((1/4-abs(x-xmodulo-1/2))>0) ...
+ ((1/4-abs(x-xmodulo+1/2))>0);
elseif (key==2)
 uex = max(0, 1-4*abs(x-xmodulo-1/2)) ...
 + max(0, 1-4*abs(x-xmodulo+1/2));
else

 uex = sin(2*pi*(x-xmodulo)) ;
 end;
if mod(it, 50)

% clf();
 subplot(1,2,1), plot(x, u, '.', x, uex, '-');
title('Discrete solution (uh)');
 subplot(1,2,2), plot(x, u-uex, '.');
 title('Error (uh-u_ex)');
 
drawnow();
end;
end;