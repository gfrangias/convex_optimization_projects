function [x,fun_val,iter,x_hist,fun_val_hist]=projected_gradient_method(f,g,h,p,x0,epsilon)
   
    x=x0;
    x_hist=x0;
    x_prev = 10^10;
    fun_val=f(x);
    fun_val_hist=fun_val;
    grad=g(x);
    hessian=h(x);
    L=max(eig(hessian));
    iter=0;
    while(norm(x-x_prev)>epsilon)
        iter=iter+1;
        x_prev=x;
        x_hist=[x_hist (x-grad/L)];
        x=p(x-grad/L);
        x_hist=[x_hist x];
        fun_val=f(x);
        fun_val_hist=[fun_val_hist fun_val];
        grad=g(x);
        hessian=h(x);
        L=max(eig(hessian));
        fprintf('iter_number = %3d norm_grad = %2.6f fun_val = %2.6f \n',...
        iter,norm(grad),fun_val);
    end
end