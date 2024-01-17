clc, close all, clear all

%%%%%%%%%%%%%%%%%%C.(a)%%%%%%%%%%%%%%%%%%%%
n=2;
m=20;
A=randn(m,n);
b=rand(m,1);
c=randn(n,1);

cvx_begin
    variable x(n)
    minimize( c'*x-sum(log(b-A*x)) )
    subject to
        min(b-A*x)>0;
cvx_end

fprintf('\n Press a key to continue...\n\n'),
pause

%%%%%%%%%%%%%%%%%%C.(b)%%%%%%%%%%%%%%%%%%%%
if n==2
    x_opt = x;
    x_1 = x_opt-0.5:0.01:x_opt+0.5;
    x_2 = x_1;
    lenx = length(x_1);

    for ii = 1:lenx
        for jj = 1:lenx
            x = [x_1(ii);x_2(jj)];
            if min(b-A*x)>0
                f(ii,jj) = c'*x-sum(log(b-A*x));
            else
                f(ii,jj) = 35;
            end
        end
    end

    figure(1)
    mesh(x_1, x_2, f'), xlabel('x1'), ylabel('x2'), zlabel('f(x1,x2)'), hold on
    clc, fprintf('\nI plot the function. Press a key to continue...\n'), pause

    figure(2)
    contour(x_1, x_2, f',10), xlabel('x1'), ylabel('x2'), title('Level sets of the function')
    clc, fprintf('\nI plot the level sets of the function. Press a key to continue...\n'), pause
end

%%%%%%%%%%%%%%%%%%C.(c)%%%%%%%%%%%%%%%%%%%%

x0=zeros(n,1); 
s=1;        % always 1
alpha=0.5;  % from 0 to 0.5
beta=0.1;   % from 0 to 1
epsilon=0.1;
f = @(x) c'*x-sum(log(b-A*x));
g = @(x) diag(c+sum(A./(b-A*x)));
constr = @(x) min(b-A*x);

fprintf('\n Backtracking Line Search ...\n\n');
[x_b,fun_val_b,iter_b,x_hist_b,fun_val_hist_b] = ...
    gradient_method_backtracking_constrained(f,g,x0,s,alpha,beta,epsilon,constr);
fprintf('\n Backtracking Line Search ^\n\n');

fprintf('\n Press a key to continue...\n\n'),
pause

%%%%%%%%%%%%%%%%%%C.(d)%%%%%%%%%%%%%%%%%%%%

h = @(x) A'*A./((b-A*x)'*(b-A*x));

fprintf('\n Newton''s Method ...\n\n');
[x_n,fun_val_n,iter_n,x_hist_n,fun_val_hist_n] = ...
    newton_method(f,g,h,x0,s,alpha,beta,epsilon,constr);
fprintf('\n Newton''s Method ^\n\n');

fprintf('\n Press a key to continue...\n\n'),
pause

%%%%%%%%%%%%%%%%%%C.(e)%%%%%%%%%%%%%%%%%%%%

k_b = 0:iter_b;
k_n = 0:iter_n;

dim_str = strcat(num2str(n),'x',num2str(m));
figure(3)
semilogy(k_b,fun_val_hist_b-cvx_optval,"b.-",k_n,fun_val_hist_n-cvx_optval,"r.-"),
xlabel('k'), ylabel('f(x_k)-p_\ast'), legend('Backtracking Line Search','Newton''s Method'),
title(['Dimensions: ',+dim_str]),
clc, fprintf('\nI plot the level sets of the function.\n')

