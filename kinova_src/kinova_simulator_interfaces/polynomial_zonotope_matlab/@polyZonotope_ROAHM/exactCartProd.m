function pZ = exactCartProd(pZ1,pZ2)
% by Bohao Zhang?? (Patrick Holmes copied and pasted)
% cartProd - Returns the cartesian product of two polyZonotope object with
% same factors
%
% Syntax:  
%    pZ = cartProd(pZ1,Z2)
%
% Inputs:
%    pZ1 - polyZonotope object
%    pZ2 - polyZonotope or contSet object
%
% Outputs:
%    pZ - polyZonotope object
%
% Example: 
%    pZ = polyZonotope(2,[1 3 1],[],[1,2,3]);
%    zono = zonotope([1,3]);
%
%    pZcart = cartProd(pZ,zono);
%
%    plot(pZcart,[1,2],'r','Filled',true,'EdgeColor','none');
%    xlim([1 8]);
%    ylim([-3 5]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/cartProd

% Author:       Niklas Kochdumper
% Written:      25-June-2018
% Last update:  05-May-2020 (MW, standardized error message)
% Last revision:---

%------------- BEGIN CODE --------------
    
    % convert other set representations to polyZonotopes (first set)
    if isa(pZ1,'double') || isa(pZ2,'double')
        pZ = cartProd(pZ1, pZ2);
        return;
    end
    
    if isa(pZ1,'zonotope') || isa(pZ2,'zonotope')
        pZ = cartProd(pZ1, pZ2);
        return;
    end
    
    if ~isa(pZ1,'polyZonotope_ROAHM') || ~isa(pZ2,'polyZonotope_ROAHM')
        error('exactCartProd only supports two polyZonotope_ROAHM')
    end

    % get dimensions
    dim1 = length(pZ1.c);
    dim2 = length(pZ2.c);
    
    % center vector
    c = [pZ1.c;pZ2.c];
    
    % generator matrix, exponent matrix and identifier vector
    if isempty(pZ1.G)
        if isempty(pZ2.G)
           G = []; 
           expMat = []; 
           id = [];
        else
           G = [zeros(dim1,size(pZ2.G,2));pZ2.G];
           expMat = pZ2.expMat;
           id = pZ2.id;
        end
    else
        if isempty(pZ2.G)
           G = [pZ1.G;zeros(dim2,size(pZ1.G,2))];
           expMat = pZ1.expMat;
           id = pZ1.id;
        else
           [id,expMat1,expMat2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.expMat,pZ2.expMat);
            
           G = blkdiag(pZ1.G,pZ2.G);
           expMat = [expMat1,expMat2];
        end
    end
    
    % matrix of independent generators
    Grest = [];
    
    if isempty(pZ1.Grest)
       if ~isempty(pZ2.Grest)
           Grest = [zeros(dim1,size(pZ2.Grest,2));pZ2.Grest];
       end
    else
       if isempty(pZ2.Grest)
           Grest = [pZ1.Grest;zeros(dim2,size(pZ1.Grest,2))];
       else
           Grest = blkdiag(pZ1.Grest,pZ2.Grest);
       end
    end       
        
    % generate new polyZonotope
    pZ = polyZonotope_ROAHM(c,G,Grest,expMat,id);
    pZ = deleteZeros(pZ);
%     pZ.initialId = pZ1.initialId;

%------------- END OF CODE --------------