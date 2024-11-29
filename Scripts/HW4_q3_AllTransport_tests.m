format("shortE");
sz = [3 5];
varTypes = {'string', 'double', 'double', 'double', 'double'};
varNames = {'IC', 'a=1,nu=1/2', 'a=-1,nu=1/2', 'a=1,nu=1.01', 'a=-1,nu=1'};
T = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames);
T(:,1) = {'S'; 'SF'; 'HF'};

A = [1 -1 1 -1];
NU = [1/2 1/2 1.01 1];
THETA = [-1 0 1];
keys = [0 1 2];

N = 10^3;
NumIts = 10^4;
h = 1/N; % spatial grid size
x = h/2 : h : 1-h/2; % discretization of the interval [0,1]

% Computes all possible combinations of
% k -> theta = -1, 0, or 1
% j -> a=1,nu=1/2 or a=-1,nu=1/2 or a=1,nu=1.01 or a=-1,nu=1
% i -> key = 0, 1, or 2
for k = 1:3
    theta = THETA(k);
    scheme = theta+1;
    for j = 1:4
        a = A(j);
        nu = NU(j);
        lambda = nu/abs(a); % CFL number
        dt = lambda*h; % temporal grid size
        for i = 1:3
            key = keys(i);
            % Initial condition u0
            if key == 1
                u = 0 + ((1/4-abs(x-1/2))>0);
            elseif key == 2
                u = max(0, 1-4*abs(x-1/2));
            else
                u = sin(2*pi*x);
            end
            
            % Compute solution using conservation form of PDE
            time = 0;
            for it = 1:NumIts
                time = time + dt;
                % computing numerical flux function
                Phi = 0.5*a*([u, u(1)] + [u(N), u]) ...
                - 0.5*nu^theta*abs(a)*([u, u(1)] - [u(N), u]);
                % updating u with forward time-step
                u = u - lambda * (Phi(2:N+1) - Phi(1:N));
            
                % computing exact solution for error comparison
                xmodulo = mod(a*time, 1); % need modulo for periodic boundaries
                if key == 1
                    uex = 0 + ((1/4-abs(x-xmodulo-1/2)) > 0) ...
                    + ((1/4-abs(x-xmodulo+1/2)) > 0);
                elseif key == 2
                    uex = max(0, 1-4*abs(x-xmodulo-1/2)) ...
                    + max(0, 1-4*abs(x-xmodulo+1/2));
                else
                    uex = sin(2*pi*(x-xmodulo));
                end
            end
            
            f = figure();
            % solution
            subplot(1,2,1), plot(x, u, '.', x, uex, '-');
            title('Discrete solution (u_h)');
            % error
            subplot(1,2,2), plot(x, u-uex, '.');
            title('Error (u_h - u_{ex})');
            
            % saving
            formatSpec = '/Figures/HW4_q3_end_%d%d%d.png';
            str = sprintf(formatSpec, scheme, i, j);
            exportgraphics(f, [pwd str]);
            close;
    
            % appending table
            T(i,j+1) = {round(norm(u-uex), 4, 'significant')};
        end
    end
    formatSpec = '/Tables/HW4_q3_table%d.xlsx';
    tbl = sprintf(formatSpec, scheme);
    writetable(T, [pwd tbl], 'Sheet', 1);
end