classdef (InferiorClasses = {?polyZonotope}) polyZonotope_ROAHM < polyZonotope
    % ROAHM version of polynomial zonotopes
    % patrick holmes 20210901
    % copied from polyZonotope.m by Niklas Kochdumper

% polyZonotope - Object Constructor for polonomial zonotopes 
%
% Syntax:  
%    object constructor:    Obj = polyZonotope(c,G,Grest)
%                           Obj = polyZonotope(c,G,Grest,expMat)
%                           Obj = polyZonotope(c,G,Grest,expMat,id)
%
% Inputs:
%    c - center of the polynomial zonotope (dimension: [nx,1])
%    G - generator matrix containing the dependent generators 
%       (dimension: [nx,N])
%    Grest - generator matrix containing the independent generators
%            (dimension: [nx,M])
%    expMat - matrix containing the exponents for the dependent generators
%             (dimension: [p,N])
%    id - vector containing the integer identifiers for the dependent
%         factors (dimension: [p,1])
%
% Outputs:
%    Obj - Polynomial Zonotope Object
%
% Example: 
%    c = [0;0];
%    G = [2 0 1;0 2 1];
%    Grest = [0;0.5];
%    expMat = [1 0 3;0 1 1];
%
%    pZ = polyZonotope(c,G,Grest,expMat);
%
%    plot(pZ,[1,2],'r','Filled',true,'EdgeColor','none');  
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope

% Author:        Niklas Kochdumper
% Written:       26-March-2018 
% Last update:   02-May-2020 (MW, add property validation, def constructor)
% Last revision: ---

%------------- BEGIN CODE --------------
    methods
        function obj = polyZonotope_ROAHM(varargin)
            obj = obj@polyZonotope(varargin{:});
        end
        
        % methods in separate files
        c = cross(a, b);
        c = plus(a, b);
        c = minus(a, b);
        c = cos(a, b);
        c = sin(a, b);
        c = times(a, b);
        c = power(a, b);
        c = exactCartProd(a, b);
        pZ = reduce(pZ,option,order,varargin);
        [c, d] = pz_to_msspoly(a, b);
        [c] = pz_to_sym(a, b);
        c = slice(a, b);
        c = remove_dependence(a, b)
    end
end

