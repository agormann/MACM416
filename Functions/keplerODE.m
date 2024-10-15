function dydt = keplerODE(t,y)

dydt = zeros(4, 1);
dydt(1) = y(3);
dydt(2) = y(4);
dydt(3) = -y(1) ./ sqrt((y(1).^2 + y(2).^2).^3);
dydt(4) = -y(2) ./ sqrt((y(1).^2 + y(2).^2).^3);

end

