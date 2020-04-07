% Script to perform excercise 2.4
clear;clc;close all;

% We define a handle for the Beale function
f = @(x) (1.5 - x(1) + x(1)*x(2))^2 + (2.25 - x(1) + x(1)*x(2)^2)^2 + (2.625 - x(1) + x(1)*x(2)^3)^2;
stepsize = 0.01;
[X,Y] = meshgrid(-3:stepsize:3);
z = f(X,Y);
niveles = [0.1, 3:5];
contour(X,Y,z, niveles)

axis equal

% example iteration history
hold on
x = [2:-0.2:0];
y = x/2;
plot(x,y,'--d')
hold off




