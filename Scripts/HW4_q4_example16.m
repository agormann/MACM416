% example16: spatial differences. 
% We compare 4 different methods for computing the second derivative.
clear;
close all;
%% Problem Specification
% x in [-pi, pi] with periodic boundary conditions, i.e.
% x_0 = -pi = pi = x_N

M = 11; % maximum power of 2
N = 10*2.^(0:M); % list of trial N values

% function to compute the derivative of as well as its analytic solution
phi = @(x) -cos(x); true = @(x) cos(x);
% phi = @(x) cos(2*x); true = @(x) -4*cos(2*x);

% cell arrays for storing results
cell_cpu = cell(M+1, 4);
cell_err = cell(M+1, 4);
cell_err_max = cell(M+1, 4);
cell_grid = cell(M+1);
%% Computations
for i = 1:M+1
    n = N(i);
    h = 2*pi/n; % mesh width
    x = h*(1:n)'; % mesh. recall -pi = pi, so [-pi,pi] ~ [0,2pi]
    cell_grid{i} = x; % storing for plotting x vs. error
    
    % method 1
    % computing the first order (centered) derivative twice
    % we replace u'(x) by (u(x_{i+1}) - u(x_{i-1})) / 2h
    tic; % timing how long the method takes to build
    firstColumn = [0 -1 zeros(1,n-2)]'; 
    firstRow = [0  1 zeros(1,n-2)];
    SD = toeplitz(firstColumn, firstRow); % forming spatial differencing matrix
    SD(1,n) = -1; % enforcing periodic boundary conditions
    SD(n,1) =  1;
    SD = SD/(2*h);
    SD2 = SD*SD;
    cell_cpu{i,1} = toc; % storing time
    
    err_SD2 = SD2*phi(x) - true(x); % computing error
    cell_err{i,1} = err_SD2; % storing result
    cell_err_max{i,1} = max(err_SD2); % storing max result

    % method 2
    % computing second order (centered) derivative
    % we replace u''(x) by (u(x_{i+1}) - 2u(x_i) + u(x_{i-1})) / h^2
    tic;
    firstColumn = [-2 1 zeros(1,n-2)];
    TD = toeplitz(firstColumn);
    TD(1,n) = 1;
    TD(n,1) = 1;
    TD = TD/h^2;
    cell_cpu{i,2} = toc;
    
    err_TD = TD*phi(x) - true(x);
    cell_err{i,2} = err_TD;
    cell_err_max{i,2} = max(err_TD);
    
    % method 3
    % spectral differentiation (explicit)
    % leveraging fourier analysis to compute the derivative; inverse DFT
    % 1. we expand the periodic grid in terms of shifted delta functions
    % which are interpolated by periodic sinc function
    % 2. then compute the derivative (twice) in this representation
    tic;
    firstColumn = [0 .5*(-1).^(1:n-1).*cot((1:n-1)*h/2)]; % derivative of sinc
    firstRow = firstColumn([1 n:-1:2]); % reversed of firstColumn
    SpecD = toeplitz(firstColumn, firstRow);
    SpecD = SpecD*SpecD;
    W = SpecD*phi(x);
    cell_cpu{i,3} = toc;

    err_SpecD = W - true(x);
    cell_err{i,3} = err_SpecD;
    cell_err_max{i,3} = max(err_SpecD);

    % method 4
    % spectral differentiation (FFT)
    % 1. compute DFT of u -> U
    % 2. take second derivative in fourier space -> W_k = (ik)^2 U_k = -k^2 U_k
    % 3. compute iDFT of W_k -> w_k
    tic;
    U = fft(phi(x)); % DFT of u(x)
    W = ifft(-[0:n/2 1-n/2:-1]'.^2 .* U); % second derivative
    cell_cpu{i,4} = toc;
    
    err_FFT = W - true(x);
    cell_err{i,4} = err_FFT;
    cell_err_max{i,4} = max(err_FFT);
end
%% Plotting

% making xticklabels as {'10', '20', '40', ...}
xx = 1:M+1;
labels = cell(1,M+1);
for i = 1:M+1
    labels{1,i} = sprintf('%d', N(i));
end

%% cpu time vs. N
f1 = figure();
legend('Location', 'best');
grid on;

xlabel('Size of Matrix (N)');
xticks(1:M+1);
xticklabels(labels);
xlim([1 M+1]);
ylabel('CPU Time (s)');

hold on;
plot(xx', [cell_cpu{:,1}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#0072BD", ...
    'DisplayName', 'SD2');
plot(xx', [cell_cpu{:,2}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#D95319", ...
    'DisplayName', 'TD');
plot(xx', [cell_cpu{:,3}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#EDB120", ...
    'DisplayName', 'Spectral');
plot(xx', [cell_cpu{:,4}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#7E2F8E", ...
    'DisplayName', 'FFT');
%% max absolute error vs. N
f2 = tiledlayout(1,2); % SD2 & TD
xlabel(f2, 'Size of Matrix (N)');
ylabel(f2, 'Max Absolute Error');

nexttile;
start = 1;
legend('Location', 'best');
grid on;
xticks(1:M+1);
xticklabels(labels);
xlim([1 M+1]);

hold on;
plot(xx(start:end), [cell_err_max{start:end,1}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#0072BD", ...
    'DisplayName', 'SD2');
plot(xx(start:end), [cell_err_max{start:end,2}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#D95319", ...
    'DisplayName', 'TD');

nexttile;
start = 6;
legend('Location', 'best');
grid on;
xticks(1:M+1);
xticklabels(labels);
xlim([start M+1]);

hold on;
plot(xx(start:end), [cell_err_max{start:end,1}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#0072BD", ...
    'DisplayName', 'SD2');
plot(xx(start:end), [cell_err_max{start:end,2}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#D95319", ...
    'DisplayName', 'TD');
%% 
f3 = tiledlayout(1,2); % spectral diff & FFT diff
xlabel(f3, 'Size of Matrix (N)');
ylabel(f3, 'Max Absolute Error');

nexttile;
start = 1;
legend('Location', 'best');
grid on;
xticks(1:M+1);
xticklabels(labels);
xlim([start M+1]);

hold on;
plot(xx(start:end), [cell_err_max{start:end,3}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#EDB120", ...
    'DisplayName', 'Spectral');
plot(xx(start:end), [cell_err_max{start:end,4}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#7E2F8E", ...
    'DisplayName', 'FFT');

nexttile;
start = 6;
legend('Location', 'best');
grid on;
xticks(1:M+1);
xticklabels(labels);
xlim([start M+1]);

hold on;
plot(xx(start:end), [cell_err_max{start:end,3}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#EDB120", ...
    'DisplayName', 'Spectral');
plot(xx(start:end), [cell_err_max{start:end,4}], '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#7E2F8E", ...
    'DisplayName', 'FFT');
%% (max absolute error x cpu time) vs. N
f4 = tiledlayout(1,2); % SD2 & TD
xlabel(f4, 'Size of Matrix (N)');
ylabel(f4, 'Max Absolute Error x CPU Time');

nexttile;
start = 1;
legend('Location', 'best');
grid on;
xticks(1:M+1);
xticklabels(labels);
xlim([start M+1]);

hold on;
plot(xx(start:end), ([cell_err_max{start:end,1}].*[cell_cpu{start:end,1}]), '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#0072BD", ...
    'DisplayName', 'SD2');
plot(xx(start:end), ([cell_err_max{start:end,2}].*[cell_cpu{start:end,2}]), '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#D95319", ...
    'DisplayName', 'TD');

nexttile;
start = 6;
legend('Location', 'best');
grid on;
xticks(1:M+1);
xticklabels(labels);
xlim([start M+1]);

hold on;
plot(xx(start:end), ([cell_err_max{start:end,1}].*[cell_cpu{start:end,1}]), '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#0072BD", ...
    'DisplayName', 'SD2');
plot(xx(start:end), ([cell_err_max{start:end,2}].*[cell_cpu{start:end,2}]), '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#D95319", ...
    'DisplayName', 'TD');
%%
f5 = tiledlayout(1,2); % spectral diff & FFT diff
xlabel(f5, 'Size of Matrix (N)');
ylabel(f5, 'Max Absolute Error x CPU Time');

nexttile;
start = 1;
legend('Location', 'best');
grid on;
xticks(1:M+1);
xticklabels(labels);
xlim([start M+1]);

hold on;
plot(xx(start:end), ([cell_err_max{start:end,3}].*[cell_cpu{start:end,3}]), '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#EDB120", ...
    'DisplayName', 'SD2');
plot(xx(start:end), ([cell_err_max{start:end,4}].*[cell_cpu{start:end,4}]), '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#7E2F8E", ...
    'DisplayName', 'TD');

nexttile;
start = 6;
legend('Location', 'best');
grid on;
xticks(1:M+1);
xticklabels(labels);
xlim([start M+1]);

hold on;
plot(xx(start:end), ([cell_err_max{start:end,3}].*[cell_cpu{start:end,3}]), '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#EDB120", ...
    'DisplayName', 'SD2');
plot(xx(start:end), ([cell_err_max{start:end,4}].*[cell_cpu{start:end,4}]), '.-', 'MarkerSize', 10, 'LineWidth', 1, 'Color', "#7E2F8E", ...
    'DisplayName', 'TD');