%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M file that computes and plots the Taylor approximations of Exercise_1 %
%                                                                        %
% A. P. Liavas, October 19, 2022                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all

% (1) First and second-order Taylor of f(x) = 1/(1+x), in R_+.
x(:,1)=[0:.1:5];  % construct the x_axis
x_0 = 5;         % Taylor point - you may change this point

% Compute function f at the points of the x_axis
f = 1./(1+x);     
% Plot the function
figure(1)
plot(x, f), xlabel('x axis'), ylabel('f(x)'), hold on, grid on
fprintf('I plot the function. Press a key to continue...'), pause

% Compute first-order Taylor approximation around x_0 
f1 = 1/(1+x_0) - (1+x_0)^(-2) * (x-x_0);                                     
% Add to the plot the 1-st order approximation
plot(x, f1, 'r'), ylabel('f(x) and 1-st order Taylor')
clc, fprintf('\nI plot the function and its 1-st order approximation. Press a key to continue...'), pause

% Compute second-order Taylor approximation around x_0
f2 = 1/(1+x_0) - (1+x_0)^(-2) * (x-x_0) + 1/2 * 2 * (1+x_0)^(-3) * (x-x_0).^2;   
% Add to the plot the 2-nd order approximation
plot(x, f2), ylabel('f(x) and 1-st and 2-nd order Taylor') 
title('Function and its 1-st and 2-nd order Taylor approximations')
legend('f(x)', '1-st order Taylor', '2-nd order Taylor')
clc, fprintf('\nI plot the function and its 1-st and 2-nd order approximations. Press a key to continue...')

clear f f1 f2
clc
