function [x,fun_val,iter,x_hist,fun_val_hist]=gradient_method_exact(f,g,x0,P,epsilon)
% Gradient method with exact line search
%
% INPUT
%=======================================
% f ......... objective function
% g ......... gradient of the objective function
% x0......... initial point
% P ......... table P of the quadratic function
% epsilon ... tolerance parameter for stopping rule
% OUTPUT
%=======================================
% x ......... optimal solution (up to a tolerance)
% of min f(x)
% fun_val ... optimal function value
% iter ...... number of iterations needed
% x_hist .... history of optimal points estimations
% fun_val ... history of function values estimations

x=x0;
x_hist=x0;
grad=g(x);
fun_val=f(x);
fun_val_hist=fun_val;
iter=0;
while(norm(grad)>epsilon)
    iter=iter+1;
    t=norm(grad)^2/(grad'*P*grad);
    x=x-t*grad;
    x_hist=[x_hist x];
    grad=g(x);
    fun_val=f(x);
    fun_val_hist=[fun_val_hist fun_val];
    fprintf('iter_number = %3d norm_grad = %2.6f funval = %2.6f \n',...
            iter, norm(grad), fun_val);
end
end