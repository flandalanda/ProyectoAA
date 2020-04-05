clear;   close all;   clc;

%% Define f with two arguments  and  so that
% it can be evaluated for matrices of values X, Y
f = @(x) (1.5 - x(1) + x(1).*x(2)).^2 + (2.25 - x(1) + x(1).*x(2).^2).^2 + (2.625 - x(1) + x(1).*x(2).^3).^2;


%% define point and trust region radius
x0    = [0;0];
delta = 1;


%% plot f in cartesian coordinates arround x0
showPlot = true;
if showPlot
	Delta = 1.1*delta;
	X = linspace(x0(1)-Delta, x0(1)+Delta, 32);
	Y = linspace(x0(2)-Delta, x0(2)+Delta, 32);
    [X1, X2] = meshgrid(X,Y);
	Z  = arrayfun(@(x1,x2) f([x1,x2]), X1, X2);
	s2 = surf(X,Y,Z);
end


%% plot quadratic model in trust region with polar coordinates arround x0
hold on
% polar coordinates arround x0
[T,R] = meshgrid(linspace(0,2*pi,64),linspace(0,delta,16));
X     = R.*cos(T) +x0(1);
Y     = R.*sin(T) +x0(2);
% the quadratic model (simple in this case)
B = apHess(f,x0);
g = apGrad(f,x0);
m  = @(x) f(x0) + g'*x' + x*B*x';
% evaluation and plot
Z  = arrayfun(@(x1,x2) m([x1,x2]), X, Y);
s1 = mesh(X,Y,Z);
view(130, 15)
hold off
