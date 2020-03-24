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
eta = 0.1;
tol = 10^-5;
deltaMax = 1.5;
msg = "El mínimo fue encontrado de forma exitosa";

iter = 0;
x = x0;
n = length(x0);
delta = deltaMax;
grad = apGrad(f,x);
hess = apHess(f,x);

while(norm(grad)>tol && iter <itmax)
    
    l1 = min(eigs(hess))
    if l1 <= 0
        hess = hess + (10^-12-1.125*l1)*eye(n);
    end
    
    dogLeg = pDogLeg(hess,grad,delta);
    x1 = x + dogLeg;
    df = f(x)-f(x1);
    dm = -grad'*dogLeg-.5*dogLeg'*hess*dogLeg;
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
iter
if iter == itmax
    msg = "Se exedió el número de iteraciones permitidas";
end

end

