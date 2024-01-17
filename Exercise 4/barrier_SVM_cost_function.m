function [result] = barrier_SVM_cost_function( w, X_augm, y, t )
%BARRIER_SVM_COST_FUNCTION Return the result of the cost function gk(w)
%
    N = size(y,2);
    i = 1;
    result = 0.5*t*(norm(w))^2;
    while(i<=N)
        result = result - log(y(:,i)*w'*X_augm-1);
        i = i+1;
    end
end

