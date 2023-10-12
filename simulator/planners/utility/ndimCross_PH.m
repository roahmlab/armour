function [v] = ndimCross_PH(Q)
% ndimCross - Computes the n-dimensional cross product
%
% Syntax:  
%    [o] = orth(V)
%
% Inputs:
%    V - matrix of column vectors
%
% Outputs:
%    o - orthogonal vector
%
% Example: 
%    Z=zonotope(rand(2,5));
%    P=polytope(Z);
%    plot(P);
%    hold on
%    plot(Z);
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: intervalhull,  vertices

% Author: Matthias Althoff
% Written: 14-September-2006 
% Last update: 22-March-2007
% Last revision: ---

%------------- BEGIN CODE --------------

% for i=1:length(Q(:,1))
%     D=Q;
%     D(i,:)=[];
%     v(i,1)=(-1)^(i+1)*det(D);
% end

v = [Q(2, 1)*Q(3, 2) - Q(3, 1)*Q(2, 2); Q(3, 1)*Q(1, 2) - Q(1, 1)*Q(3, 2); Q(1, 1)*Q(2, 2) - Q(2, 1)*Q(1, 2)];

% if sum((v-v2).^2) > 1e-16
%     error('Patrick 3D cross is wrong.')
% end

%------------- END OF CODE --------------