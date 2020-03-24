function [pC] = pCauchy(B,g, delta)
% In: B     ... (symmetric matrix) approximate of the hessian of f at xk
%     g     ... (vector) gradient of f in xk
%     delta ... trust region radius
%
%Out: pC    ... The Cauchy point

normg = norm(g);
termCuad = g'*B*g;

p = -delta/normg*g;
disp(p)

tau = 1

if termCuad > 0
    tau = min(normg^3/(delta*termCuad),1)
end 

pC = tau * p

end

