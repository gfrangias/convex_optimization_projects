function [result] = point_is_feasible(w_int, X_augm, y)
%POINT_IS_FEASIBLE This function checks if w_int is in the domain of phi(w)
%   INPUTS:
%   w_int
%   X_augm
%   y
%   OUTPUT:
%   result -> boolean
    N = size(y,2);
    i=1;
    result = 1;
    while(i<=N)
        if(y(:,i)*w_int'*X_augm(:,i)<=1)
            result = 0;
            break;
        end
        i = i+1;
    end
end

