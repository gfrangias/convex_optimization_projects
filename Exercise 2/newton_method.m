function [x,fun_val,iter,x_hist,fun_val_hist]=newton_method(f,g,h,x0,s,alpha,beta,epsilon,constr)
% Gradient method with backtracking stepsize rule with constraint
%
% INPUT
%=======================================
% f ......... objective function
% g ......... gradient of the objective function
% h ......... hessian of the objective function
% x0......... initial point
% s ......... initial choice of stepsize
% alpha ..... tolerance parameter for the stepsize selection
% beta ...... the constant in which the stepsize is multiplied
% at each backtracking step (0<beta<1)
% epsilon ... tolerance parameter for stopping rule
% constr .... constraint of the objective function
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
hessian=h(x);
fun_val_hist=fun_val;
newton_s=-inv(hessian)*grad;    %Newton step
dec=grad'*inv(hessian)*grad;    %Newton decrement
iter=0;
while (dec/2>epsilon)
    iter=iter+1;
    t=s;
    while (constr(x-t*grad)<=0)
        t=beta*t;
        fprintf('Out of domain! iter_number = %3d new step = %2.6f\n',...
        iter,t);        
    end
    while (fun_val-f(x-t*grad)<alpha*t*norm(grad)^2)
        t=beta*t;
    end
    x=x+t*newton_s;
    x_hist=[x_hist x];
    fun_val=f(x);
    fun_val_hist=[fun_val_hist fun_val];
    grad=g(x);
    hessian=h(x);
    newton_s=-inv(hessian)*grad;
    dec=grad'*inv(hessian)*grad;
    
    fprintf('iter_number = %3d decrement/2 = %2.6f fun_val = %2.6f \n',...
    iter,0.5*dec,fun_val);
end