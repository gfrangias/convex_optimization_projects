clc, close all; 
clear;

%%%%%%%%%%%%%%%%%%5.(a)%%%%%%%%%%%%%%%%%%%%
n=500;
p=100;
x0=10*randn(n,1);
A=randn(p,n);
b=randn(p,1);

cvx_begin
    variable x(n)
    minimize( norm(x-x0) )
    subject to
        A*x-b==0;
cvx_end

opt = x0-A'*inv(A*A')*(A*x0-b);
opt_val = norm(x0-opt);
dist = sqrt(sum((x-x0).^2));
dist2 = sqrt(abs((x0-opt)'*A'*inv(A*A')*(A*x0-b)));
clc,fprintf('Distance: %2.6f My distance: %2.6f\n', dist, dist2), pause

%%%%%%%%%%%%%%%%%%5.(b)%%%%%%%%%%%%%%%%%%%%
clc,clear;
n=2;
p=1;
A = randn(p,n);
%rank(A)
b=randn(p,1);

% Creating positive definite matrix P
B = rand(n,n);

[U,S,V] = svd(B);

K = 2;
l_min = rand(1,1);
l_max = K*l_min;

z = l_min + (l_max-l_min)*rand(n-2,1);
eig_P = [l_min;l_max;z];
Lambda = diag(eig_P);

q = randn(n,1);
P = U * Lambda * U';

% CVX solution
cvx_begin
    variable x(n)
    minimize( 0.5*x'*P*x+q'*x )
    subject to
        A*x-b==0;
cvx_end
pause
x_opt_cvx = x;

% Solving the linear equations system
factor1 = [P,A';A,zeros(p)];
product = [-q;b];
factor2 = linsolve(factor1,product);
primal_optimal = factor2(1:n);
dual_optimal = factor2(n+1:n+p);
opt_val = 0.5*primal_optimal'*P*primal_optimal+q'*primal_optimal;

clc,fprintf('CVX''s opt_val: %2.6f My opt_val: %2.6f\n', cvx_optval, opt_val)

% Projected Gradient Method
f = @(x) 0.5*x'*P*x+q'*x;
g= @(x) P*x+q;
h= @(x) P;
p = @(x) x-A'*inv(A*A')*(A*x-b);
epsilon=0.01;
x0=randn(n,1);

[x,fun_val,iter,x_hist,fun_val_hist]=projected_gradient_method(f,g,h,p,x0,epsilon);

if(n==2)
    max_dist = max(abs(primal_optimal(1,:)-x0(1)),abs(primal_optimal(2,:)-x0(2)));
    x1 = (primal_optimal(1,:)-max_dist-1:0.01:primal_optimal(1,:)+max_dist)+1;
    x2 = (primal_optimal(2,:)-max_dist-1:0.01:primal_optimal(2,:)+max_dist)+1;
    for ii=1:length(x1)
        for jj=1:length(x2)
            x3 = [x1(ii);x2(jj)];
            f2(ii,jj) = 0.5*x3'*P*x3+q'*x3;
        end
    end
    y=(-A(1)*x1 + b)/A(2);
    projections = x_hist(:,3:2:end);
    figure(1)
    hold on;
    axis equal
    contour(x1, x2, f2'), xlabel('x1'), ylabel('x2'), title('Level sets of the function')
    plot(x_hist(1,:),x_hist(2,:),"b.-")
    plot(projections(1,:),projections(2,:),"r*")
    plot(x1,y);
end




