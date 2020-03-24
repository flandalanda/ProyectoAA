function [p] = pDogLeg( B, g, delta )
% In :  B     ... an s.p.d. matrix that approximates the hessian of  f  in xk
%       g     ... (vector) gradient of  f  in  xk
%       delta ... trust region radius
%
% Out:  p     ... The  dogleg  point
%
% Algorithm taken from class notes uploaded to Comunidad ITAM.

% We calculate the constants we will need throughout the algorithm
normg = norm(g);
const = delta/normg;

% We get the minimizer of m along the steepest descent direction
alpha = normg^2/(g'*B*g);

% We verify the minimizer remains in the trust region, else we take the
% largest step possible in that direction
if alpha >= const
    p = -const*g;
else
    % We calculate the Newton step which is the unconstrained minimizer of
    % the quadratic model.
    pB = -B\g;
    
    % If the full Newton step remains in the trust region, we take it
    if norm(pB) <= delta
        p = pB;
    else
        % We look for a point in the line segment connecting pU and pB
        % We seek the positive solution to ||pU+(tau-1)(pB-pU)||^2 = delta^2
        pU = -0.99*alpha*g;
        aux = pB-pU;

        % Upon simplification the problem is expressed as a quadratic
        % scalar function on tau, which can be solved explicitly
        a = dot(aux,aux);
        b = 2*dot(pU, aux);
        c = dot(pU,pU)-delta^2;

        % We obtain the positive root of the quadratic function
        tau = 0.5*(-b + sqrt(b^2-4*a*c))/a;
        
        % We define our step
        p = pU + (tau-1)*aux;
    end
end

end

