function [O] = make_orientation(theta, axis)

switch axis
    case 1
        O = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
    case 2
        O = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
    case 3
        O = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    otherwise
        error('incorrect axis spec.');
end

end

