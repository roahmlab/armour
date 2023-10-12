function [pZ] = plus(pZ1, pZ2)
	% NOTE: this implements the equivalent of @polyZonotope/exactPlus,
	% which is different from @polyZonotope/plus

	% patrick holmes 20210901
	% overloading plus for polyZonotope_ROAHM
	% adapted from @polyZonotope/exactPlus by Niklas Kochdumper

% exactPlus - compute the addition of two sets while preserving the
%             dependencies between the two sets
%
% Syntax:  
%    pZ = exactPlus(pZ1,pZ2)
%
% Inputs:
%    pZ1 - polyZonotope object
%    pZ2 - polyZonotope object
%
% Outputs:
%    pZ - polyZonotope object after exact addition
%
% Example: 
%    pZ1 = polyZonotope([0;0],[2 1 2;0 2 2],[],[1 0 3;0 1 1]);
%    pZ2 = [1 2;-1 1]*pZ1;
%   
%    pZ = pZ1 + pZ2;
%    pZ_ = exactPlus(pZ1,pZ2);
%
%    figure
%    subplot(1,2,1);
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none','Splits',10);
%    title('Minkowski Sum');
%    subplot(1,2,2);
%    plot(pZ_,[1,2],'b','Filled',true,'EdgeColor','none');
%    title('Exact Addition');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes, zonotope/plus

% Author:       Niklas Kochdumper
% Written:      26-March-2018 
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------
    
if isa(pZ1, 'polyZonotope_ROAHM') && isa(pZ2, 'polyZonotope_ROAHM')
    if isempty(pZ1.G) && isempty(pZ2.G)
        pZ = polyZonotope_ROAHM(pZ1.c + pZ2.c, [], [pZ1.Grest,pZ2.Grest]);
    else
        % bring the exponent matrices to a common representation
        [id,expMat1,expMat2] = mergeExpMatrix(pZ1.id,pZ2.id,pZ1.expMat,pZ2.expMat);

        % add up all generators that belong to identical exponents
        [ExpNew,Gnew] = removeRedundantExponents([expMat1,expMat2],[pZ1.G,pZ2.G]);

        % assemble the properties of the resulting polynomial zonotope
        pZ = polyZonotope_ROAHM(pZ1.c + pZ2.c, Gnew, [pZ1.Grest,pZ2.Grest], ExpNew, id);
    end  
elseif isnumeric(pZ1) && isa(pZ2, 'polyZonotope_ROAHM')
    pZ = polyZonotope_ROAHM(pZ1 + pZ2.c, pZ2.G, pZ2.Grest, pZ2.expMat, pZ2.id);
elseif isnumeric(pZ2) && isa(pZ1, 'polyZonotope_ROAHM')
    pZ = polyZonotope_ROAHM(pZ2 + pZ1.c, pZ1.G, pZ1.Grest, pZ1.expMat, pZ1.id);
elseif isnumeric(pZ1) && isnumeric(pZ2)
    pZ = pZ1 + pZ2;
end
end

