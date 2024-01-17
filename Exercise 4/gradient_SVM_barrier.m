function [result] = gradient_SVM_barrier(w, X_augm, y, t)
%GRADIENT_SVM_BARRIER Return the gradient of the cost function gk(w)
%   
    N = size(y,2);
    i = 1;
    result = t*w;
    while(i<=N)
        result = result + (y(:,i)/(1-y(:,i)*X_augm(:,i)'*w))*X_augm(:,i);
        i = i+1;
    end
end

