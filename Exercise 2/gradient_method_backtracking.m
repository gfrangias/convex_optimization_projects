function [x,fun_val,iter,x_hist,fun_val_hist]=gradient_method_backtracking(f,g,x0,s,alpha,...
beta,epsilon)
% Gradient method with backtracking stepsize rule
%
% INPUT
%=======================================
% f ......... objective function
% g ......... gradient of the objective function
% x0......... initial point
% s ......... initial choice of stepsize
% alpha ..... tolerance parameter for the stepsize selection
% beta ...... the constant in which the stepsize is multiplied
% at each backtracking step (0<beta<1)
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
while (norm(grad)>epsilon)
iter=iter+1;
t=s;
while (fun_val-f(x-t*grad)<alpha*t*norm(grad)^2)
t=beta*t;
end
x=x-t*grad;
x_hist=[x_hist x];
fun_val=f(x);
fun_val_hist=[fun_val_hist fun_val];
grad=g(x);
fprintf('iter_number = %3d norm_grad = %2.6f fun_val = %2.6f \n',...
iter,norm(grad),fun_val);
end