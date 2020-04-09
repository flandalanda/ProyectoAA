% Script to perform excercise 2.4
clear;clc;close all;

% We define a handle for the Beale function
f = @(x,y) (1.5 - x + x.*y).^2 + (2.25 - x + x.*y.^2).^2 + (2.625 - x + x.*y.^3).^2;

% Definimos el punto inicial, el radio de la Region de Confianza y las
% variables de entrada
x0 = [-3; -3];
xstar = [3;0.5];
delta = 0.5;
itmax = 1000;

%Definimos los conjuntos de nivel
stepsize =  0.01;
[X,Y] = meshgrid(-5:stepsize:5);
Z = f(X,Y);
niveles = [0.1, 1:100];
contour(X,Y,Z,niveles)

axis equal
hold on

%Ejecuta ejercicio 2.3
Ej2_3

%Plot Punto Cauchy

Xc = table1(:,5);
    
plot(Xc(1),Xc(2),'--d')
hold on    
    
%Plot Punto DogLeg

Xdl = table2(:,5);
    
plot(Xdl(1),Xdl(2),'--d')
