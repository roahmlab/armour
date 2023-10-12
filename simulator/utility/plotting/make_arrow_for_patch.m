function [F,V] = make_arrow_for_patch(p,varargin)   
% [faces,vertices] = make_arrow_for_patch(p_end)
% [faces,vertices] = make_arrow_for_patch(p_start,p_end)
% [faces,vertices] = make_arrow_for_patch(p_start,'property1',value1,'property2',value2,...)
% [faces,vertices] = make_arrow_for_patch(p_start,p_end,'property1',value1,'property2',value2,...)
%
% This function generates faces and vertices that can be used with patch
% to plot a 3-D arrow from the origin to p_end, or else from a point
% p_start to the point p_end; all such points should be passed in as
% 3-by-1 vectors (in \R^3).
%
% This function can also take in three different property/value pairs:
%   'shaft_width'   the diameter of the cylindrical arrow shaft
%   'head_length'   the length of the conical arrowhead
%   'head_width'    the diameter of the conical arrowhead
%
% EXAMPLE USAGE:
%   % make arrow from origin to (1,1,1)
%   [F,V] = make_arrow_for_patch(ones(3,1))
%
%   % plot the arrow using patch
%   figure(1) ; clf ; grid on ; hold on ; axis equal ; view(3)
%   patch('faces',F,'vertices',V,'facealpha',0.1)
%
% See also: plot_arrow
 
    if mod(length(varargin),2) == 1
        p_start = p ;
        p_end = varargin{1} ;
        
        if length(varargin) > 1
            varargin = varargin(2:end) ;
        else
            varargin = {} ;
        end
    else
        p_start = zeros(3,1) ;
        p_end = p ;
    end
    
    % set default values based on how far apart p1 and p2 are
    d = norm(p_end - p_start) ;
    shaft_width = d/20 ;
    head_length = d/5 ;
    head_width = d/10 ;
    
    % iterate through varargin
    for idx = 1:2:length(varargin)
        switch lower(varargin{idx})
            case 'shaft_width'
                shaft_width = varargin{idx+1} ;
            case 'head_length'
                head_length = varargin{idx+1} ;
            case 'head_width'
                head_width = varargin{idx+1} ;
            otherwise
                error('Invalid keyword!')
        end
    end

    % make the shaft and shift it up
    shaft_length = d - head_length ;
    [F_shaft,V_shaft] = make_cylinder_for_patch(shaft_width/2,shaft_length) ;
    V_shaft(:,3) = V_shaft(:,3) + shaft_length/2 ;
    
    % make the head and shift it up
    [F_head,V_head] = make_cone_for_patch(head_length,head_width/2) ;
    V_head(:,3) = V_head(:,3) + shaft_length ;
    
    % stack all the vertices
    N_shaft = size(V_shaft,1) ;
    F_head = F_head + N_shaft ;
    F = [F_shaft ; F_head] ;
    V = [V_shaft ; V_head] ;
    
    % get the rotation matrix for the direction between p1 and p2
    n0 = [0;0;1] ;
    n = (p_end - p_start)./d ;
    v = cross(n0,n) ;
    if ~all(v == 0)
        v_hat = skew(v) ;
        s = norm(v) ;
        c = dot(n0,n) ;
        R = eye(3) + v_hat + v_hat^2 .* ((1-c)./(s^2)) ;
    else
        R = eye(3) ;
    end
    
    % rotate and translate the vertices
    V = (R*V')' + repmat(p_start(:)',size(V,1),1) ;
end