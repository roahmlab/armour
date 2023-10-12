function pZ = times(factor1,factor2)
    % patrick holmes 20211103
    % adapted from @polyZonotope/mtimes by Niklas Kochdumper
    % multiplication of a PolyZonotope_ROAHM with a polyZonotope_ROAHM
    % *both have to be 1-dimensional*
    % need to keep track of generator coefficients.

% times - Overloaded '*' operator for the multiplication of a polyZonotope 
%   with a polyZonotope
%
% Syntax:  
%    pZ = mtimes(pZ1,pZ2)
%
% Inputs:
%    factor1 - polyZonotope_ROAHM object or scalar
%    factor2 - polyZonotope_ROAHM object or scalar
%
% Outputs:
%    pZ - polyZonotope_ROAHM after multiplication of scalars
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus, zonotope/mtimes

% Author:       Niklas Kochdumper
% Written:      25-June-2018 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% if isa(factor1,'PolyZonotope_ROAHM')
%     pZ1=factor1;
% else
%     error('first arg must be matPolyZonotope_ROAHM');
% end
% 
% if isnumeric(factor2)
%     pZ2 = factor2;
%     c_new = pZ1.C*pZ2;
%     G_new = [];
%     Grest_new = [];
%     expMat_new = [];
%     id = pZ1.id;
%     if ~isempty(pZ1.G)
%         Gc = pagemtimes(pZ1.G, pZ2);
%         G_new = [G_new, Gc(:, :)];
%         expMat_new = [expMat_new, pZ1.expMat];
%     end
%     if ~isempty(pZ1.Grest)
%         Gcrest = pagemtimes(pZ1.Grest, pZ2);
%         Grest_new = [Grest_new, Gcrest(:, :)];
%     end
        
if isa(factor1, 'polyZonotope_ROAHM') && isa(factor2, 'polyZonotope_ROAHM')
    pZ1 = factor1;
    pZ2 = factor2;
    if pZ1.dim ~= 1 || pZ2.dim ~=1
        error('both polyZonotope_ROAHM objects must have dimension 1 to use this times function');
    end
    % bring the exponent matrices to a common representation
    [id,expMat1,expMat2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.expMat,pZ2.expMat);

    % % add up all generators that belong to identical exponents
    % [ExpNew,Gnew] = removeRedundantExponents([expMat1,expMat2],[pZ1.G,pZ2.G]);

    % multiply 
    G_new = [];
    Grest_new = [];
    expMat_new = [];
    % get new center
    c_new = pZ1.c*pZ2.c;
    % deal with dependent gens
    if ~isempty(pZ2.G)
        cg = pZ1.c*pZ2.G;
        G_new = [G_new, cg];
        expMat_new = [expMat_new, expMat2];
    end
    if ~isempty(pZ1.G)
        Gc = pZ1.G*pZ2.c;
        G_new = [G_new, Gc];
        expMat_new = [expMat_new, expMat1];
    end
    if ~isempty(pZ1.G) && ~isempty(pZ2.G)
        GG = (pZ1.G'*pZ2.G)';
        G_new = [G_new, GG(:)'];
        for i = 1:size(pZ1.G, 2)
           expMat_new = [expMat_new, expMat1(:, i) + expMat2];
        end
    end
    % deal with independent gens
    if ~isempty(pZ2.Grest)
        cgrest = pZ1.c*pZ2.Grest;
        Grest_new = [Grest_new, cgrest];
    end
    if ~isempty(pZ1.Grest)
        Gcrest = pZ1.Grest*pZ2.c;
        Grest_new = [Grest_new, Gcrest(:, :)];
    end
    if ~isempty(pZ1.Grest) && ~isempty(pZ2.Grest)
        GGrest = (pZ1.Grest'*pZ2.Grest)';
        Grest_new = [Grest_new, GGrest(:)'];
    end
    % deal with mixed gens (dependent times indenpendent, treat as independent)
    if ~isempty(pZ1.G) && ~isempty(pZ2.Grest)
        GGrest = (pZ1.G'*pZ2.Grest)';
        Grest_new = [Grest_new, GGrest(:)'];
    end
    if ~isempty(pZ1.Grest) && ~isempty(pZ2.G)
        GGrest = (pZ1.Grest'*pZ2.G)';
        Grest_new = [Grest_new, GGrest(:)'];
    end
elseif isa(factor1, 'polyZonotope_ROAHM') && isnumeric(factor2)
    pZ = factor1;
    c_new = factor2*pZ.c;
    G_new = factor2*pZ.G;
    Grest_new = factor2*pZ.Grest;
    expMat_new = pZ.expMat;
    id = pZ.id;
elseif isnumeric(factor1) && isa(factor2, 'polyZonotope_ROAHM')
    pZ = factor2;
    c_new = factor1*pZ.c;
    G_new = factor1*pZ.G;
    Grest_new = factor1*pZ.Grest;
    expMat_new = pZ.expMat;
    id = pZ.id;
else
    error('first & second arg must be numeric or polyZonotope_ROAHM');
end

% create output
if isempty(G_new) && isempty(Grest_new)
    pZ = polyZonotope_ROAHM(c_new);
elseif isempty(G_new)
    pZ = polyZonotope_ROAHM(c_new, G_new, Grest_new);
else
    pZ = polyZonotope_ROAHM(c_new, G_new, Grest_new, expMat_new, id);
end


