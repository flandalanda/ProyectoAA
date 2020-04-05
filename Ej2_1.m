clear;   close all;   clc;

%% Define the Beale function to be evaluated at a vector x
f = @(x) (1.5 - x(1) + x(1).*x(2)).^2 + (2.25 - x(1) + x(1).*x(2).^2).^2 + (2.625 - x(1) + x(1).*x(2).^3).^2;


%% Define point and trust region radius
x0    = [0;0];
delta = 1;


%% Plot f in cartesian coordinates around x0
showPlot = true;
if showPlot
	Delta = 1.1*delta;
	X = linspace(x0(1)-Delta, x0(1)+Delta, 32);
	Y = linspace(x0(2)-Delta, x0(2)+Delta, 32);
    [X1, X2] = meshgrid(X,Y);
	Z  = arrayfun(@(x1,x2) f([x1,x2]), X1, X2);
	s2 = surf(X,Y,Z);
end


%% Plot quadratic model in the trust region with polar coordinates around x0
hold on
% Polar coordinates arround x0
[T,R] = meshgrid(linspace(0,2*pi,64),linspace(0,delta,16));
X     = R.*cos(T) +x0(1);
Y     = R.*sin(T) +x0(2);
% The quadratic model
B = apHess(f,x0);
g = apGrad(f,x0);
m  = @(x) f(x0) + g'*x' + x*B*x';
% Evaluation and plot
Z  = arrayfun(@(x1,x2) m([x1,x2]), X, Y);
s1 = mesh(X,Y,Z);
view(130, 15)
hold off
