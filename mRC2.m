function [x, msg] = mRC2( f, x0, itmax )
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

% Initializing initial aproximations and iteration count
iter = 0;
x = x0;
n = length(x0);
delta = deltaMax;
grad = apGrad(f,x);
hess = apHess(f,x);

% The loop stops when the maximum number of iterations is reached or when
% we are close enough to a stationary point
while(norm(grad)>tol && iter <itmax)
    
    % We obtain the minimum eigenvector of hessian matrix
    l1 = min(eigs(hess))
    % If the minimum eigenvector is less than or equal to 0, the matrix is
    % not positve definite and we add a multiple of the identity matrix to
    % make it s.p.d.
    if l1 <= 0
        hess = hess + (10^-12-1.125*l1)*eye(n);
    end
    
    % DogLeg point is calculated
    dogLeg = pDogLeg(hess,grad,delta);
    
    % Actual reduction is calculated
    x1 = x + dogLeg;
    df = f(x)-f(x1);
    
    % Predicted reduction is calculated as mk(0)-mk(dogLeg)
    dm = -grad'*dogLeg-.5*dogLeg'*hess*dogLeg;
    
    % We calculate the reduction quotient
    rho = df/dm;
    
    % If the model fits the function poorly, we reduce the trust region
    if rho < .25
       delta = 0.25*delta; 
    % If the model is a good fit and we took the largest step possible, we
    % extend the trust region
    elseif rho > .75 && norm(dogLeg) == delta
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
    
    iter = iter +1;
end
% Output to console the number of iterations taken
iter
% Verification that loop finished becuase maximum number of iterations was
% reached or we got close enough to a stationary point
if iter == itmax
    msg = "Se exedió el número de iteraciones permitidas";
end

end

