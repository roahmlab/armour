function [C,V] = make_circle(r,n,c)
% C = make_circle(r,n,p)
% [F,V] = make_circle(r,n,p)
%
% Make a circle of radius r, approximated as an n-sided polygon, centered
% at c \in \R^2. By default, r = 1, n = 100, and c = zeros(2,1). The output
% C is a 2-by-n array.
%
% The optional second output is a list of vertices so that you can plot the
% circle with patch('faces',F,'vertices',V,other_patch_args)
%
% Authors: Shreyas Kousik
% Created: who knows!
% Updated: 7 Jan 2021

    if nargin < 3
        c = zeros(2,1) ;
    end

    if nargin < 2
        n = 100 ;
    end
    
    if nargin < 1
        r = 1 ;
    end
    
    % make angle vector
    a_vec = linspace(0,2*pi,n) ;
    
    % make C
    C = [r*cos(a_vec) ; r*sin(a_vec) ] + repmat(c(:),1,n) ;
    
    if nargout > 1
        V = C' ;
        C = [1:n, 1] ;
    end
end