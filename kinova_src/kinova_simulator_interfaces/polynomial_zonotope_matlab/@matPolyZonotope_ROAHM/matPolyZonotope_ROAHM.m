classdef (InferiorClasses = {?polyZonotope_ROAHM, ?polyZonotope}) matPolyZonotope_ROAHM
    % ROAHM matrix version of polynomial zonotopes
    % patrick holmes 20210825
    % based off polyZonotope.m by Niklas Kochdumper
    %
    % note that G(:, :, 1) is the first matrix generator, G(:, :, 2) is the
    % second, and so on.
    
properties (SetAccess = protected, GetAccess = public)
    C (:,:,1) {mustBeNumeric,mustBeFinite} = [];
    G (:,:,:) {mustBeNumeric,mustBeFinite} = [];
    Grest (:,:,:) {mustBeNumeric,mustBeFinite} = [];
    expMat (:,:) {mustBeInteger} = [];
    id (:,1) {mustBeInteger} = [];
end
    
methods
    
    function Obj = matPolyZonotope_ROAHM(C,G,Grest,varargin)
        
        if nargin == 1 
            % Copy Constructor
            if isa(C,'matPolyZonotope_ROAHM')
                Obj = C;
                
            % Single point
            elseif ismatrix(C)
                Obj.C = C;
            end
        
        % No exponent matrix provided
        elseif nargin == 2 || nargin == 3
            
            Obj.C = C;
            Obj.G = G;
            if nargin == 3
                Obj.Grest = Grest;
            end
            
            % construct exponent matrix under the assumption
            % that all generators are independent
            Obj.expMat = eye(size(G,3));
            Obj.id = (1:size(G,3))';
            
        
        % Exponent matrix as user input
        elseif nargin == 4 || nargin == 5
            
            expMat = varargin{1};
            
            % check correctness of user input
            if ~all(all(floor(expMat) == expMat)) || ~all(all(expMat >= 0)) || ...
               (~isempty(G) && size(expMat,2) ~= size(G,3))
                error('Invalid exponent matrix!');
            end
            
            % remove redundant exponents
%             if ~isempty(expMat)
%                 [expMat,G] = removeRedundantExponents(expMat,G);
%             end
            
            % assign properties
            Obj.C = C;
            Obj.G = G;
            Obj.Grest = Grest;
            Obj.expMat = expMat;
            
            % vector of integer identifiers
            if nargin == 5
                
                id = varargin{2};
                
                % check for correctness
                if length(id) ~= size(expMat,1)
                   error('Invalid vector of identifiers!'); 
                end
                
                Obj.id = id;
                
            else
                Obj.id = (1:size(expMat,1))';
            end
        

        % Else if not enough or too many inputs are passed    
        elseif nargin > 5
            disp('This class needs less input values.');
            Obj=[];
        end
        
    end
    
    % extra functions
    pZ = reduce(pZ,option,order,varargin);
    pZ = polyZonotope_ROAHM(mpz);

end

% methods (Static = true)
%     Z = generateRandom(varargin) % generate random polyZonotope
% end


end

