% Script to perform excercise 2.2
clear;   close all;   clc;

%% Define the Beale function to be evaluated at a vector x
f = @(x) (1.5 - x(1) + x(1).*x(2)).^2 + (2.25 - x(1) + x(1).*x(2).^2).^2 + (2.625 - x(1) + x(1).*x(2).^3).^2;

% The quadratic model
x0=[0;0];
B = apHess(f,x0);
g = apGrad(f,x0);
m  = @(x,y) f(x0) + g(1).*x + g(2).*y + B(1,1)*x.^2 + B(2,2)*y.^ + 2.*B(1,2).*x.*y;

%Descent directions
pdog=pDogLeg(B,g,4);
pCau=pCauchy(B,g,4);
pNewton=inv(B)*g;

%level sets
stepsize =  0.01;  % if we make smaller level sets, we get more detail 
[X,Y] = meshgrid(-5:stepsize:5);
z = m(X,Y);
niveles = [0.1, 1:100];
contour(X,Y,z, niveles)

axis equal

%trust region boundry 
hold on
ang=0:0.01:2*pi; 
xp=4*cos(ang);
yp=4*sin(ang);
plot(xp,yp);

% Plotting descent directions, dogleg is red, Cauchy is blue and Newton is
% green
quiver( 0, 0, pdog(1),pdog(2),'r','LineWidth',2);
quiver(0,0, pCau(1),pCau(2),'b', 'LineWidth',2);
quiver(0,0, pNewton(1),pNewton(2),'g', 'LineWidth',2);
