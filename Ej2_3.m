% Script to perform excercise 2.3
clear;clc;

% We define a handle for the Beale function
f = @(x) (1.5 - x(1) + x(1)*x(2))^2 + (2.25 - x(1) + x(1)*x(2)^2)^2 + (2.625 - x(1) + x(1)*x(2)^3)^2;

%We define x0 as the vector (-3,-3)' and verify the Hessian (or at
%least its approximation) is not positive definite as well as the distance
%between x0 and the known optimum
x0 = [-3;-3];
xstar = [3;0.5];
B = apHess(f,x0);

%We define 100 as the maximum number iterations
itermax = 100;

assert(norm(x0-xstar)>1.5, 'The starting point is too close to the known optimum')
assert(min(eigs(B)) < 0,'The Hessian approximation is positive definite at this point')

%We find the solution using the trust region algorithm that employs the
%Cauchy point
[sol1, msg1, table1] = mRC1(f,x0,itermax)

%We find the solution using the trust region algorithm that utilizes the
%dogleg method
[sol2, msg2, table2] = mRC2(f,x0,itermax)