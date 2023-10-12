function [pZ]=reduce(pZ,option,order,varargin)
% updated by Patrick Holmes 20211103
% if we need to make changes do so here...
% updated 20220407 to work with matPolyZonotopes
%
% reduce - Reduces the order of a polynomial zonotope
%
% Syntax:  
%    [pZ]=reduce(pZ,option,order)
%
% Inputs:
%    pZ - polyZonotope object
%    option - reduction algorithm (see zonotope/reduce)
%    order - order of reduced polynomial zonotope
%
% Outputs:
%    pZ - reduced polynomial zonotope
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 1 1],[0.1,-0.4;0.2,0.3],[1 0 3;0 1 1]);
%    pZred = reduce(pZ,'girard',2);
%
%    hold on
%    plot(pZred,[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/reduce

% Author:       Niklas Kochdumper
% Written:      23-March-2018 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    % extract dimensions
    N = size(pZ.C,1);
    N2 = size(pZ.C,2);
    P = size(pZ.G,3);
    Q = size(pZ.Grest,3);

    % number of generators that stay unreduced (N generators are added again
    % after reduction)
    K = N*N2*order - N*N2;

    % check if it is necessary to reduce the order
    if P + Q > N*order && K >= 0

        % concatenate all generators
%         G = [pZ.G,pZ.Grest];
        G = cat(3, pZ.G, pZ.Grest);

        % half the generator length for exponents that are all even
        temp = prod(ones(size(pZ.expMat))-mod(pZ.expMat,2),1);
        ind = find(temp == 1);
        G(:, :, ind) = 0.5 * G(:, :, ind);

        % calculate the length of the generator vectors with a special metric
        len = sum(sum(G.^2,1), 2);
        len = len(:);

        % determine the smallest generators (= generators that are removed)
        [~,ind] = sort(len,'descend');
        ind = ind(K+1:end);

        % split the indizes into the ones for dependent and independent
        % generators
        indDep = ind(ind <= P);
        indInd = ind(ind > P);
        indInd = indInd - P * ones(size(indInd));

        % construct a zonotope from the generators that are removed
        Grem = pZ.G(:,:,indDep);
        Erem = pZ.expMat(:,indDep);

        GrestRem = pZ.Grest(:,:,indInd);

        mpZtemp = matPolyZonotope_ROAHM(zeros(N,N2),Grem,GrestRem,Erem);
        pZtemp = polyZonotope_ROAHM(mpZtemp);

        zono = zonotope(pZtemp);    % zonotope over-approximation

        % reduce the constructed zontope with the reduction techniques for
        % linear zonotopes
        zonoRed = reduce(zono,option,1,varargin{:});

        % remove the generators that got reduced from the generator matrices
        pZ.G(:,:,indDep) = [];
        pZ.expMat(:,indDep) = [];
        pZ.Grest(:,:,indInd) = [];

        % add the reduced generators as new independent generators 
        pZ.C = pZ.C + reshape(center(zonoRed), N, N2);
        g_tmp = generators(zonoRed);
        
        pZ.Grest = cat(3, pZ.Grest, reshape(g_tmp, N, N2, size(g_tmp, 2)));

    end

    % remove all exponent vector dimensions that have no entries
    temp = sum(pZ.expMat,2);
    ind = find(temp > 0);
    pZ.expMat = pZ.expMat(ind,:);
    pZ.id = pZ.id(ind);
    
%     if pZ.dim == 1
%         pZ.Grest = sum(abs(pZ.Grest));
%     end

end


%------------- END OF CODE --------------