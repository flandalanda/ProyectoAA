% Script to perform excercise 2.4
clear;clc;close all;

% We define a handle for the Beale function
f = @(x) (1.5 - x(1) + x(1)*x(2))^2 + (2.25 - x(1) + x(1)*x(2)^2)^2 + (2.625 - x(1) + x(1)*x(2)^3)^2;

% Definimos el punto inicial, el radio de la Region de Confianza y las
% variables de entrada
x0 = [-3; -3];
xstar = [3;0.5];
delta = 0.5;
itmax = 1000;

%Definimos los conjuntos de nivel
stepsize =  0.01;
[X,Y] = meshgrid(0:stepsize:5);
Z = f(X,Y);
niveles = [0.1, 1:10];
contour(X,Y,Z, niveles)

axis equal
hold on
