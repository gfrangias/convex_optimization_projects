% Exercise 6.(b)

clear, clc, close all

m=3;
n=2;
A = randi([-20 20],m,n);
x = randi([-20 20],n,1);
b=A*x;
x1=(x(1)-3:.1:x(1)+3);
x2=(x(2)-3:.1:x(2)+3);
len_x = length(x1);

for ii=1:len_x
    for jj=1:len_x
        x3 = [x1(ii);x2(jj)];
        v = A*x3 - b;
        f(ii,jj) = v(1)^2 + v(2)^2 + v(3)^2;
    end
end

figure(7)
mesh(x1, x2, f'), xlabel('x1'), ylabel('x2'), zlabel('f(x1,x2)'), hold on
clc, fprintf('\nI plot the function f. Press a key to continue...'), pause

figure(8)
contour(x1, x2, f'), xlabel('x1'), ylabel('x2'), title('Level sets of the function')
clc, fprintf('\nI plot the level sets of the function. Press a key to continue...'), pause

clear b
clear f

e = randi([-10,10],m,1);
%e = [0;0;0];
b=A*x+10^(-2)*e;

for ii=1:len_x
    for jj=1:len_x
        x3 = [x1(ii);x2(jj)];
        v = A*x3 - b;
        f(ii,jj) = v(1)^2 + v(2)^2 + v(3)^2;
    end
end

figure(9)
mesh(x1, x2, f'), xlabel('x1'), ylabel('x2'), zlabel('f(x1,x2)'), hold on
clc, fprintf('\nI plot the function f. Press a key to continue...'), pause

figure(10)
contour(x1, x2, f'), xlabel('x1'), ylabel('x2'), title('Level sets of the function')
clc, fprintf('\nI plot the level sets of the function. Press a key to continue...'), pause

