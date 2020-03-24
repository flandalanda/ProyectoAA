function [p] = pDogLeg( B, g, delta )
% In :  B     ... an s.p.d. matrix that approximates the hessian of  f  in xk
%       g     ... (vector) gradient of  f  in  xk
%       delta ... trust region radius
%
% Out:  p     ... The  dogleg  point

% We calculate the constants we will need throughout the algorithm
normg = norm(g);
const = delta/normg;

% We get the minimizer of m along the steepest descent direction
alpha = normg^2/(g'*B*g);

if alpha >= const
    p = -const*g;
else
    % We calculate the Newton step
    pB = -B\g;
    if norm(pB) <= delta
        p = pB;
    else
        % We seek the positive solution to ||pU+(tau-1)(pB-pU)||^2 = delta^2
        pU = -0.99*alpha*g;
        aux = pB-pU;

        a = dot(aux,aux);
        b = 2*dot(pU, aux);
        c = dot(pU,pU)-delta^2;

        tau = 0.5*(-b + sqrt(b^2-4*a*c))/a;
        
        p = pU + (tau-1)*aux;
    end
end

end

