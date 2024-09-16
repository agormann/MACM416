%%Example code 2
% This is based on several codes from the book 'Spectral methods in Matlab'
% by L.N. Trefethen. 
 
f=@(x)exp(sin(x)); fprime = @(x)cos(x).*f(x);

%%

%%%% Part 1: Fourth order differencing
% For various N, set up grid in [-pi,pi] and function u(x):
  Nvec = 2.^(3:12);
 figure(1);
 clf, subplot('position',[.1 .4 .8 .5])
 tic;
  for N = Nvec
    h = 2*pi/N; x = -pi + (1:N)'*h; % (A) WHAT IS HAPPENING HERE?
    u=f(x);uprime=fprime(x); % (B) WHAT'S BEING DONE HERE?
    
    % Construct sparse fourth-order differentiation matrix:
    e = ones(N,1);
    D =   sparse(1:N,[2:N 1],2*e/3,N,N)...
        - sparse(1:N,[3:N 1 2],e/12,N,N);
    D = (D-D')/h;

    % Plot max(abs(D*u-uprime)):
    error = norm(D*u-uprime,inf); %(C) WHAT DO  YOU THINK IS HAPPENING HERE?
    loglog(N,error,'.','markersize',15), hold on
  end
  Time1=toc;
  grid on, xlabel N, ylabel error
  title('Convergence of fourth-order finite differences')
  semilogy(Nvec,Nvec.^(-4),'--') 
  text(105,5e-8,'N^{-4}','fontsize',18)


  %%Part 2: Spectral Differentiation, version 1
   
% For various N (even), set up grid as before:
 figure(2), clf, subplot('position',[.1 .4 .8 .5])
 tic;
  for N = 2:2:100;
    h = 2*pi/N;
    x = -pi + (1:N)'*h;
    u = f(x); uprime = fprime(x); 

    % Construct spectral differentiation matrix:
    column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)];
    D = toeplitz(column,column([1 N:-1:2]));

    % Plot max(abs(D*u-uprime)):
    error = norm(D*u-uprime,inf);
    loglog(N,error,'.','markersize',15), hold on
  end
  Time2=toc;
  grid on, xlabel N, ylabel error
  title('Convergence of spectral differentiation')

  %% Part 3: spectral differentiation, part 2

  % For various N (even), set up grid as before:
  figure(3), clf, subplot('position',[.1 .4 .8 .5])
  tic
  for N = 2:2:100;
    h = 2*pi/N;
    x = -pi + (1:N)'*h;
    u = f(x); uprime = fprime(x);

    %fast spectral differentiation. HOW IS THIS DIFFERENT THAN THE PREVIOUS
    %ONE?
     v_hat = fft(u);
  w_hat = 1i*[0:N/2-1 0 -N/2+1:-1]' .* v_hat;
  w = real(ifft(w_hat));

    

    % Plot max(abs(w-uprime)):
    error = norm(w-uprime,inf);
    loglog(N,error,'.','markersize',15), hold on
  end
  Time3= toc;
  grid on, xlabel N, ylabel error
  title('Convergence of spectral differentiation')
