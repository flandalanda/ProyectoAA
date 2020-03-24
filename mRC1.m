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
eta = 0.1;
tol = 10^-5;
deltaMax = 1.5;
msg = "El mínimo fue encontrado de forma exitosa";

iter = 0;
x = x0;
delta = deltaMax;
grad = apGrad(f,x);
hess = apHess(f,x);

while(norm(grad)>tol && iter <itmax)
    
    pC = pCauchy(hess,grad,delta);
    x1 = x + pC;
    df = f(x)-f(x1);
    dm = -grad'*pC-.5*pC'*hess*pC;
    rho = df/dm;
    
    if rho < .25
       delta = 0.25*delta; 
    elseif rho > .75 && norm(x1) == delta
        delta = min(2*delta, deltaMax);
    end
    
    if rho > eta
        x = x1;
    end
    
    grad = apGrad(f,x);
    hess = apHess(f,x);
    
    iter = iter +1;
end
if iter == itmax
    msg = "Se exedió el número de iteraciones permitidas";
end

end

