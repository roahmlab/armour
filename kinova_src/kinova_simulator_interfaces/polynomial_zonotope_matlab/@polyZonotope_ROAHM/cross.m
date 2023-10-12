function [c] = cross(a, b)
	% patrick holmes 20210825
	% `cross` function for polynomial zonotopes.
    if isnumeric(a) 
        A = [0, -a(3), a(2); ...
            a(3), 0, -a(1); ...
            -a(2), a(1), 0];
    else
        c = a.c;
        g = a.G;
        grest = a.Grest;
        Z = [c, g];
        G = [];
        Grest = [];
        for j = 1:size(Z, 2)
            z = Z(:, j);
            % https://en.wikipedia.org/wiki/Cross_product#Conversion_to_matrix_multiplication
            M = [0, -z(3), z(2); ...
                 z(3), 0, -z(1); ...
                 -z(2), z(1), 0];
            if j == 1
                C = M;
            else
                G = cat(3, G, M);
            end
        end
        
        for j = 1:size(grest, 2)
            z = grest(:, j);
            M = [0, -z(3), z(2); ...
                 z(3), 0, -z(1); ...
                 -z(2), z(1), 0];
            Grest = cat(3, Grest, M);
        end
        A = matPolyZonotope_ROAHM(C, G, Grest, a.expMat, a.id);
    end
	c = A*b;
end

