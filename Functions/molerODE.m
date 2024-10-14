function dydt = molerODE(t, y, year, rate)

% constants
d = 8.64;
mu_1 = 4.95*10^2;
mu_2 = 4.95*10^-2;
v_s = 0.12;
v_d = 1.23;
w = 10^-3;
k_1 = 2.19*10^-4;
k_2 = 6.12*10^-5;
k_3 = 0.997148;
k_4 = 6.79*10^-2;

% variables
h_s = 1/k_3 * (y(2) - sqrt(y(2).^2 - k_3*y(4).*(2*y(2)-y(4))));
c_s = 1/2 * (y(4)-h_s);
p_s = 1/c_s * k_4*h_s.^2;

% source
f = pchiptx(year, rate, t);

% pdes
dydt = zeros(5,1);
dydt(1) = 1/d*(p_s-y(1)) + f./mu_1;
dydt(2) = 1/v_s * (w*(y(3)-y(2)) - k_1 - mu_2/d*(p_s-y(1)));
dydt(3) = 1/v_d * (k_1 - w*(y(3)-y(2)));
dydt(4) = 1/v_s * (w*(y(5)-y(4)) - k_2);
dydt(5) = 1/v_d * (k_2 - w*(y(5)-y(4)));

end

