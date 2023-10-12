function taulist = SkewInverseDynamics(thetalist, dthetalist, ddthetalist, ...
                                   g, Ftip, Mlist, Glist, Slist)
% *** CHAPTER 8: DYNAMICS OF OPEN CHAINS ***
% Takes thetalist: n-vector of joint variables,
%       dthetalist: n-vector of joint rates,
%       ddthetalist: n-vector of joint accelerations,
%       g: Gravity vector g,
%       Ftip: Spatial force applied by the end-effector expressed in frame 
%             {n+1},
%       Mlist: List of link frames {i} relative to {i-1} at the home 
%              position,
%       Glist: Spatial inertia matrices Gi of the links,
%       Slist: Screw axes Si of the joints in a space frame, in the format
%              of a matrix with the screw axes as the columns.
% Returns taulist: The n-vector of required joint forces/torques.
% This function uses forward-backward Newton-Euler iterations to solve the 
% equation:
% taulist = Mlist(thetalist) * ddthetalist + c(thetalist, dthetalist) ...
%           + g(thetalist) + Jtr(thetalist) * Ftip

n = size(thetalist, 1);
Mi = eye(4);
Ai = zeros(6, n);
AdTi = zeros(6, 6, n + 1);
Vi = zeros(6, n + 1);
Vdi = zeros(6, n + 1);
Vdi(4: 6, 1) = -g;
AdTi(:, :, n + 1) = Adjoint(TransInv(Mlist(:, :, n + 1)));
Fi = Ftip;
taulist = zeros(n, 1);
for i=1: n    
    Mi = Mi * Mlist(:, :, i);
    Ai(:, i) = Adjoint(TransInv(Mi)) * Slist(:, i);    
    AdTi(:, :, i) = Adjoint(MatrixExp6(VecTose3(Ai(:, i) ...
                    * -thetalist(i))) * TransInv(Mlist(:, :, i)));    
    Vi(:, i + 1) = AdTi(:, :, i) * Vi(:, i) + Ai(:, i) * dthetalist(i);
    Vdi(:, i + 1) = AdTi(:, :, i) * Vdi(:, i) ...
                    + Ai(:, i) * ddthetalist(i) ...
                    + ad(Vi(:, i + 1)) * Ai(:, i) * dthetalist(i);    
end
for i = n: -1: 1
    
    w = crossop(Vi(1:3, i+1));
    I = Glist(1:3,1:3,i);
    m = Glist(4,4,i);
    
    Ci = [w*I + I*w 0*w; 0*w m*w];
    
    Fi = AdTi(:, :, i + 1)' * Fi ...
       + Glist(:, :, i) * Vdi(:, i + 1) ...
       + Ci * Vi(:, i + 1);
%        - ad(Vi(:, i + 1))' * (Glist(:, :, i) * Vi(:, i + 1));
     
    taulist(i) = Fi' * Ai(:, i);
end
end


function W = crossop(w)

W = [0   -w(3)   w(2);
    w(3)   0    -w(1);
   -w(2)  w(1)     0];
end