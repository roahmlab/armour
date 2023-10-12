function pZ = mtimes(factor1,factor2)
    % patrick holmes 20210825
    % adapted from @polyZonotope/mtimes by Niklas Kochdumper
    % multiplication of a matPolyZonotope_ROAHM with a polyZonotope_ROAHM
    % need to keep track of generator coefficients.

% mtimes - Overloaded '*' operator for the multiplication of a matrix or an
%          interval matrix with a polyZonotope
%
% Syntax:  
%    pZ = mtimes(matrix,pZ)
%
% Inputs:
%    matrix - numerical or interval matrix
%    pZ - polyZonotope object 
%
% Outputs:
%    pZ - polyZonotope_ROAHM after multiplication of a matrix with a polyZonotope
%
% Example: 
%    pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0;0],[1 0 3;0 1 1]);
%    matrix = [1 2;-1 3];
%    intMatrix = interval([0.9 1.9;-1.1 2.9],[1.1 2.1;-0.9 3.1]);
%       
%    pZres = matrix*pZ;
%    pZresInt = intMatrix*pZ;
%   
%    figure
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');
%    figure
%    plot(pZres,[1,2],'b','Filled',true,'EdgeColor','none');
%    figure
%    plot(pZresInt,[1,2],'g','Filled',true,'EdgeColor','none');
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

% find a matPolyZonotope_ROAHM object
if isa(factor1,'matPolyZonotope_ROAHM')
    pZ1=factor1;
else
    error('first arg must be matPolyZonotope_ROAHM');
end

if isnumeric(factor2)
    pZ2 = factor2;
    c_new = pZ1.C*pZ2;
    G_new = [];
    Grest_new = [];
    expMat_new = [];
    id = pZ1.id;
    if ~isempty(pZ1.G)
        Gc = pagemtimes(pZ1.G, pZ2);
        G_new = [G_new, Gc(:, :)];
        expMat_new = [expMat_new, pZ1.expMat];
    end
    if ~isempty(pZ1.Grest)
        Gcrest = pagemtimes(pZ1.Grest, pZ2);
        Grest_new = [Grest_new, Gcrest(:, :)];
    end
    
    % create output
    if isempty(G_new) && isempty(Grest_new)
        pZ = polyZonotope_ROAHM(c_new);
    elseif isempty(G_new)
        pZ = polyZonotope_ROAHM(c_new, G_new, Grest_new);
    else
        pZ = polyZonotope_ROAHM(c_new, G_new, Grest_new, expMat_new, id);
    end
        
elseif isa(factor2, 'polyZonotope_ROAHM')
    pZ2 = factor2;
    % bring the exponent matrices to a common representation
    [id,expMat1,expMat2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.expMat,pZ2.expMat);

    % % add up all generators that belong to identical exponents
    % [ExpNew,Gnew] = removeRedundantExponents([expMat1,expMat2],[pZ1.G,pZ2.G]);

    % multiply
    % (note, this is the same process as the multiplication implemented for
    % "rotatotopes". we're multiplying a matPolyZonotope pZ1 with center C and
    % generator matrix G, times a polyZonotope pZ2 with center c and generators
    % g. So we have to keep track of C*c (no coefficients), C*g (pZ2's
    % coefficients), G*c (pZ1's coefficients) and G*g (pZ1 and pZ2's
    % coefficients). We also have to carry out the multiplication for
    % "independent" generators in Grest.
    G_new = [];
    Grest_new = [];
    expMat_new = [];
    % get new center
    c_new = pZ1.C*pZ2.c;
    % deal with dependent gens
    if ~isempty(pZ2.G)
        Cg = pZ1.C*pZ2.G;
        G_new = [G_new, Cg];
        expMat_new = [expMat_new, expMat2];
    end
    if ~isempty(pZ1.G)
        Gc = pagemtimes(pZ1.G, pZ2.c);
        G_new = [G_new, Gc(:, :)];
        expMat_new = [expMat_new, expMat1];
    end
    if ~isempty(pZ1.G) && ~isempty(pZ2.G)
        Gg = pagemtimes(pZ1.G, pZ2.G);
        G_new = [G_new, Gg(:, :)];
        for i = 1:size(pZ1.G, 3)
           expMat_new = [expMat_new, expMat1(:, i) + expMat2];
        end
    end
    % deal with independent gens
    if ~isempty(pZ2.Grest)
        Cgrest = pZ1.C*pZ2.Grest;
        Grest_new = [Grest_new, Cgrest];
    end
    if ~isempty(pZ1.Grest)
        Gcrest = pagemtimes(pZ1.Grest, pZ2.c);
        Grest_new = [Grest_new, Gcrest(:, :)];
    end
    if ~isempty(pZ1.Grest) && ~isempty(pZ2.Grest)
        Ggrest = pagemtimes(pZ1.Grest, pZ2.Grest);
        Grest_new = [Grest_new, Ggrest(:, :)];
    end
    % deal with mixed gens (dependent times indenpendent, treat as independent)
    if ~isempty(pZ1.G) && ~isempty(pZ2.Grest)
        Ggrest = pagemtimes(pZ1.G, pZ2.Grest);
        Grest_new = [Grest_new, Ggrest(:, :)];
    end
    if ~isempty(pZ1.Grest) && ~isempty(pZ2.G)
        Ggrest = pagemtimes(pZ1.Grest, pZ2.G);
        Grest_new = [Grest_new, Ggrest(:, :)];
    end
    
    % create output
    if isempty(G_new) && isempty(Grest_new)
        pZ = polyZonotope_ROAHM(c_new);
    elseif isempty(G_new)
        pZ = polyZonotope_ROAHM(c_new, G_new, Grest_new);
    else
        pZ = polyZonotope_ROAHM(c_new, G_new, Grest_new, expMat_new, id);
    end
    
elseif isa(factor2, 'matPolyZonotope_ROAHM')
    pZ2 = factor2;
    % bring the exponent matrices to a common representation
    [id,expMat1,expMat2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.expMat,pZ2.expMat);

    % % add up all generators that belong to identical exponents
    % [ExpNew,Gnew] = removeRedundantExponents([expMat1,expMat2],[pZ1.G,pZ2.G]);

    % multiply
    % (note, this is the same process as the multiplication implemented for
    % "rotatotopes". we're multiplying a matPolyZonotope pZ1 with center C and
    % generator matrix G, times a matPolyZonotope pZ2 with center C and
    % generator matrices G. So we have to keep track of C*C (no coefficients), C*G (pZ2's
    % coefficients), G*C (pZ1's coefficients) and G*G (pZ1 and pZ2's
    % coefficients). We also have to carry out the multiplication for
    % "independent" generators in Grest.
    G_new = [];
    Grest_new = [];
    expMat_new = [];
    % get new center
    c_new = pZ1.C*pZ2.C;
    % deal with dependent gens
    if ~isempty(pZ2.G)
        Cg = pagemtimes(pZ1.C, pZ2.G);
        G_new = cat(3, G_new, Cg);
        expMat_new = [expMat_new, expMat2];
    end
    if ~isempty(pZ1.G)
        Gc = pagemtimes(pZ1.G, repmat(pZ2.C, 1, 1, size(pZ1.G, 3)));
        G_new = cat(3, G_new, Gc);
        expMat_new = [expMat_new, expMat1];
    end
    if ~isempty(pZ1.G) && ~isempty(pZ2.G)
        for i = 1:size(pZ1.G, 3)
            Gg = pagemtimes(pZ1.G(:, :, i), pZ2.G);
            G_new = cat(3, G_new, Gg);
            expMat_new = [expMat_new, expMat1(:, i) + expMat2];
        end
    end
    % deal with independent gens
    if ~isempty(pZ2.Grest)
        Cgrest = pagemtimes(pZ1.C,pZ2.Grest);
        Grest_new = cat(3, Grest_new, Cgrest);
    end
    if ~isempty(pZ1.Grest)
        Gcrest = pagemtimes(pZ1.Grest, repmat(pZ2.C, 1, 1, size(pZ1.Grest, 3)));
        Grest_new = cat(3, Grest_new, Gcrest);
    end
    if ~isempty(pZ1.Grest) && ~isempty(pZ2.Grest)
        for i = 1:size(pZ1.Grest, 3)
            Ggrest = pagemtimes(pZ1.Grest(:, :, i), pZ2.Grest);
            Grest_new = cat(3, Grest_new, Ggrest);
        end
    end
    % deal with mixed gens (dependent times indenpendent, treat as independent)
    if ~isempty(pZ1.G) && ~isempty(pZ2.Grest)
        for i = 1:size(pZ1.G, 3)
            Ggrest = pagemtimes(pZ1.G(:, :, i), pZ2.Grest);
            Grest_new = cat(3, Grest_new, Ggrest);
        end
    end
    if ~isempty(pZ1.Grest) && ~isempty(pZ2.G)
        for i = 1:size(pZ1.Grest, 3)
            Ggrest = pagemtimes(pZ1.Grest(:, :, i), pZ2.G);
            Grest_new = cat(3, Grest_new, Ggrest);
        end
    end
    
    % create output
    if isempty(G_new) && isempty(Grest_new)
        pZ = matPolyZonotope_ROAHM(c_new);
    elseif isempty(G_new)
        pZ = matPolyZonotope_ROAHM(c_new, G_new, Grest_new);
    else
        pZ = matPolyZonotope_ROAHM(c_new, G_new, Grest_new, expMat_new, id);
    end
else
    error('second arg must be numeric or polyZonotope_ROAHM or matPolyZonotope_ROAHM');
end


%-----------