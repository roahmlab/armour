function [PA, Pb, C] = polytope_PH(Z, options)
% polytope - Converts a zonotope from a G- to a H-representation
%
% This function is implemented based on Theorem 7 of
%
% Althoff, M.; Stursberg, O. & Buss, M. Computing Reachable Sets of Hybrid 
% Systems Using a Combination of Zonotopes and Polytopes Nonlinear 
% Analysis: Hybrid Systems, 2010, 4, 233-249
%
% Syntax:  
%    [P] = polytope(Z)
%
% Inputs:
%    Z - zonotope object
%
% Outputs:
%    P - polytope object
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

% Author:       Matthias Althoff
% Written:      30-September-2008
% Last update:  26-February-2009
%               05-July-2010
%               20-October-2010
%               03-December-2010
%               02-September-2011 (delete aligned added)
%               12-February-2012
%               13-October-2014
%               20-March-2015
%               02-June-2015
%               11-June-2015
%               31-August-2015
%               12-August-2016
%               30-September-2016 (one-dimensional case added)
% Last revision:---

%------------- BEGIN CODE --------------

% if nargin == 1
%     options.polytopeType = 'mpt';
% elseif nargin == 2
%     options = varargin{1};
% end

%obtain number of generators, dimensions
c = Z(:, 1);
G = Z(:, 2:end);
% G(:, ~any(G)) = [];

% patrick edit 12/03/2019
% reduce small generators for numerical reasons:
% compute metric of generators
h = vecnorm(G);
% sort generators according to metric
[h_sort, indices] = sort(h, 'descend');
threshold = 0.001;
first_reduce_idx = find(h_sort < threshold, 1, 'first');
if ~isempty(first_reduce_idx)
    Gunred = G(:, indices(1:first_reduce_idx-1));
    Gred = G(:, indices(first_reduce_idx:end));
    % box remaining generators
    d=sum(abs(Gred),2);
    %build box Gbox from interval hull vector d
    Gbox=diag(d);
    G = [Gunred, Gbox];
end
    
[dim,nrOfGenerators]=size(G);

if dim > 2
    %get number of possible facets
%     comb=combinator(nrOfGenerators,dim-1,'c');

    if ~exist('options', 'var') || nrOfGenerators > options.maxcombs
        comb=combinator(nrOfGenerators,dim-1,'c');
    else
        comb = options.combs{nrOfGenerators};
    end

    %build C matrices
    C=[];
%     for i=1:length(comb(:,1))
%         indices=comb(i,:);
%         Q=G(:,indices);
%         v=ndimCross(Q);
% %         v = ndimCross_PH(Q);
%         C(end+1,:)=v'/norm(v);
%     end
    
    Q = [G(:, comb(:, 1)); G(:, comb(:, 2))];
    if dim == 3
        C = [Q(2, :).*Q(6, :) - Q(3, :).*Q(5, :); -(Q(1, :).*Q(6, :) - Q(3, :).*Q(4, :)); Q(1, :).*Q(5, :) - Q(2, :).*Q(4, :)];
    else
        error('Dimension not supported.');
    end
    C = (C./sqrt(sum(C.^2)))';

    %remove NaN rows due to rank deficiency
    index = find(sum(isnan(C),2));
    if ~isempty(index)
        C(index,:) = [];
    end
else
    C = G;
    %... get perpendicular vector
    C = [-C(2, :); C(1, :)];
    C = (C./sqrt(sum(C.^2)))';
end

%build d vector
%determine delta d

% deltaD=zeros(length(C(:,1)),1);
% for iGen=1:nrOfGenerators
%     deltaD=deltaD+abs(C*G(:,iGen));
% end

deltaD = sum(abs((C*G)'))';

%compute dPos, dNeg
d = C*c;
% dPos=C*c+deltaD;
% dNeg=-C*c+deltaD;

%convert to mpt or ppl Polytope
% if isfield(options,'polytopeType') && strcmp(options.polytopeType,'ppl')
%     P=pplPolytope([C;-C],[dPos;dNeg]);
% else
%     P=mptPolytope([C;-C],[dPos;dNeg]);
% end

PA = [C; -C];
% Pb = [dPos; dNeg];
Pb = [d+deltaD; -d+deltaD];

%------------- END OF CODE --------------
