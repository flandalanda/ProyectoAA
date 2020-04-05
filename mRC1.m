function [x, msg] = mRC1( f, x0, itmax )
% Trust region method using the  Cauchy  point
%
% In :  f     ... (handle) function to be optimized
%       x0    ... (vector) initial point
%       itmax ... (natural number) upper bound for number of iteration
%
% Out:  x   ... (vector) last approximation of a stationary point
%       msg ... (string) message that says whether (or not) a minimum was
%       found

% We define global constraints
eta = 0.1;
tol = 10^-5;
deltaMax = 1.5;
msg = "El mínimo fue encontrado de forma exitosa";

% We initialize an empty matrix that stores the iteration, the distance
% from the current solution to the know optimum of the Beale function, the
% norm of the gradient at the current iteration and the value of the
% objective function at the current iteration in columns 1, 2, 3, and 4
% respectively

resTable(1:8,1:4)=0;

% Initializing initial aproximations and iteration count
iter = 0;
x = x0;
delta = deltaMax;
grad = apGrad(f,x);
hess = apHess(f,x);

% Create first entry in the results table
resTable(mod(iter,8)+1,:) = [iter, norm(x-[3;0.5]),norm(grad),f(x)];

% The loop stops when the maximum number of iterations is reached or when
% we are close enough to a stationary point
while(norm(grad)>tol && iter <itmax)
    
    % The Cauchy point is calculated
    pC = pCauchy(hess,grad,delta);
    
    % Actual reduction is calculated
    x1 = x + pC;
    df = f(x)-f(x1);
    % Predicted reduction is calculated as mk(0)-mk(pC)
    dm = -grad'*pC-.5*pC'*hess*pC;
    % We calculate the reduction quotient
    rho = df/dm;
    
    % If the model fits the function poorly, we reduce the trust region
    if rho < .25
       delta = 0.25*delta; 
    % If the model is a good fit and we took the largest step possible, we
    % extend the trust region
    elseif rho > .75 && norm(pC) == delta
        delta = min(2*delta, deltaMax);
    end
    
    % If none of the previous conditions were met, the trust region remains
    % the same
    
    % If the the model is a good fit for the function, we take the step
    if rho > eta
        x = x1;
        
        % If the step is taken, aproximations are recalculated
        grad = apGrad(f,x);
        hess = apHess(f,x);
    end
    
    % Increase the number of iterations performed
    iter = iter +1;
    
    % We update the latest entry in the results table
    resTable(mod(iter,8)+1,:) = [iter, norm(x-[3;0.5]),norm(grad),f(x)];
    
    
end
% Output to console the number of iterations taken
iter
% Verification that loop finished becuase maximum number of iterations was
% reached or we got close enough to a stationary point
if iter == itmax
    msg = "Se exedió el número de iteraciones permitidas";
end

% Output to console the results table
resTable

end

