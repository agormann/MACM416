function [x,res] = iterSolve(A,b,x0,tol,maxit,alg)
%ITERSOLVE Summary of this function goes here
%   Detailed explanation goes here

n = length(b);
x = x0;

% jacobi
if alg == 0         
    for k = 1:maxit
        for i = 1:n
            x(i) = 1/A(i,i) * (b(i) - A(i,:)*x0 + A(i,i)*x0(i));
        end
        if norm(x - x0) <= tol
            res = norm(x-x0);
            disp("Tolerance reached.");
            break;
        elseif k == maxit
            res = norm(x-x0);
            disp("Max iterations reached.");
            break;
        end
        x0 = x;
    end
% gauss-siedel
elseif alg == 1     
    for k = 1:maxit
        for i = 1:n
            x(i) = 1/A(i,i) * (b(i) - A(i,:)*x + A(i,i)*x(i));
        end
        if norm(x - x0) <= tol
            res = norm(x-x0);
            disp("Tolerance reached.");
            break;
        elseif k == maxit
            res = norm(x-x0);
            disp("Max iterations reached.");
            break;
        end
        x0 = x;
    end
end

end

