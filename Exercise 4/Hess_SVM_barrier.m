function [result] = Hess_SVM_barrier(w, X_augm, y, t)
%GRADIENT_SVM_BARRIER Return the hessian of the cost function gk(w)
%   
    N = size(y,2);
    i = 1;
    result = t*eye(3);
    while(i<=N)
        result = result + (1/(y(:,i)*X_augm(:,i)'*w-1)^2)*X_augm(:,i)*X_augm(:,i)';
        i = i+1;
    end
end

