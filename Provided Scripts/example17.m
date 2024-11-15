% Example 17: the 2D FFT and derivatives in 2D 
clc, clf, clear all
% Grid and initial data:EVEN GRID!! USE N even.
  N = 24; h = 2*pi/N;
x = -pi + (1:N)'*h; y = x';% what is the domain on which we work?
  [xx,yy] = meshgrid(x,y); %generating 2-D grid
 
  %%vv = exp(-40*((xx-.4).^2 + yy.^2)); %function to test
  
  vv = cos(3*xx)+sin(2*yy);% simple function to test
  figure(1), surf(xx,yy,vv), shading("flat"), view(0,90),title('u(x,y)')
  vv_truexx= -9*cos(3*xx); vv_trueyy= -4*sin(2*yy); % true derivative

  % setting up derivatives. Is there a more efficient way?
    uxx = zeros(N,N); uyy = zeros(N,N);
    ii = 1:N;
    for i = 1:N                % 2nd derivs wrt x in each row
      v = vv(i,:); V =v;
      U = (fft(V));
      W1 = (ifft(1i*[0:N/2-1 0 1-N/2:-1].*U)); % diff wrt theta
      W2 = (ifft(-[0:N/2 1-N/2:-1].^2.*U));    % diff^2 wrt theta
      uxx(i,ii) = W2(ii);
    end
    for j = 1:N                % 2nd derivs wrt y in each column
      v = vv(:,j); V =v; 
      U = (fft(V));
      W1 = (ifft(1i*[0:N/2-1 0 1-N/2:-1]'.*U));% diff wrt theta   
      W2 = (ifft(-[0:N/2 1-N/2:-1]'.^2.*U));   % diff^2 wrt theta
      uyy(ii,j) = W2(ii);
    end
    
    %% Plotting. We'll interpolate onto a finer grid
      [xxx,yyy] = meshgrid(-pi:1/48:pi,-pi:1/48:pi);
     
      figure(2)
     
 subplot(2,2,1)
      WX=uxx;vvv = interp2(xx,yy,WX,xxx,yyy,'cubic');
       surf(xxx,yyy,vvv)
      title(['uxx']),
      shading interp, view(0,90),colorbar
      %error in x derivative;
 subplot(2,2,2)
      WX=uxx-vv_truexx;vvv = interp2(xx,yy,WX,xxx,yyy,'cubic');
       surf(xxx,yyy,vvv)
      %title('Error in uxx')
      shading interp, view(0,90),colorbar
      % y derivative
       subplot(2,2,3)
      WY=uyy;vvv = interp2(xx,yy,WY,xxx,yyy,'cubic');
      surf(xxx,yyy,vvv)
      %title('uyy')
      shading interp, view(0,90),colorbar
      % error in y derivative
       subplot(2,2,4)
      WY=uyy-vv_trueyy;
       vvv = interp2(xx,yy,WY,xxx,yyy,'cubic');
      surf(xxx,yyy,vvv)
     shading interp, colorbar
    view(0,90)
   title('x')
