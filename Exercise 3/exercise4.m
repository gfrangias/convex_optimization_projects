clc, close all; 
clear;

%%%%%%%%%%%%%%%%%%4.(b)%%%%%%%%%%%%%%%%%%%%
n=2;
a=5*randn(n,1);
b=5*randn;
x0=5*randn(n,1);

cvx_begin
    variable x(n)
    minimize(0.5*(x'*x))
    subject to
        a'*x==b;
cvx_end

f = @(x) 0.5*norm(x)^2;
g= @(x) x;
h= @(x) eye(n);
p = @(x) x-a.*((a'*x-b)/(norm(a)^2));
epsilon=0.01;

[x,fun_val,iter,x_hist,fun_val_hist]=projected_gradient_method(f,g,h,p,x0,epsilon);

if(n==2)
    x_H=-10:0.1:10;
    y=(-a(1)*x_H + b)/a(2);
    figure(1)
    hold on
    axis equal
    plot(x_H,y);
    plot(x_hist(1,:),x_hist(2,:))
    grid on
end
