function [R, R_t] = get_pz_rotations_from_q(q, rotation_axis, deg)
% assuming q is a (1D) poly zonotope of angles
% use cos and sin function to convert to overapproximative
% rotation matrices.

    if ~exist('deg', 'var')
        deg = 6; % use degree 6 taylor expansion by default
    end

    cos_q = cos(q, deg);
    sin_q = sin(q, deg);
    cos_sin_q = exactCartProd(cos_q, sin_q);
    % use axis-angle formula to get rotation matrix:
    % https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    e = rotation_axis./norm(rotation_axis);
    U = [0 -e(3) e(2);...
         e(3) 0 -e(1);...
         -e(2) e(1) 0];

    % create rotation matrices from cos and sin dimensions
    cq = cos_sin_q.c(1);
    sq = cos_sin_q.c(2);
    C = eye(3) + sq*U + (1-cq)*U^2;
    C_t = C';

    Grest = [];
    Grest_t = [];
    for j = 1:size(cos_sin_q.Grest, 2)
        cq = cos_sin_q.Grest(1, j);
        sq = cos_sin_q.Grest(2, j);
        Grest_tmp = sq*U + -cq*U^2;
        Grest = cat(3, Grest, Grest_tmp);
        Grest_t = cat(3, Grest_t, Grest_tmp');
    end

    G = [];
    G_t = [];
    for j = 1:size(cos_sin_q.G, 2)
        cq = cos_sin_q.G(1, j);
        sq = cos_sin_q.G(2, j);
        G_tmp = sq*U + -cq*U^2;
        G = cat(3, G, G_tmp);
        G_t = cat(3, G_t, G_tmp');
    end

    if isempty(G) && isempty(Grest)
        R = matPolyZonotope_ROAHM(C);
        R_t = matPolyZonotope_ROAHM(C_t);
    elseif isempty(G) && ~isempty(Grest)
        R = matPolyZonotope_ROAHM(C, G, Grest);
        R_t = matPolyZonotope_ROAHM(C_t, G_t, Grest_t);
    elseif ~isempty(G)
        R = matPolyZonotope_ROAHM(C, G, Grest, cos_sin_q.expMat, cos_sin_q.id);
        R_t = matPolyZonotope_ROAHM(C_t, G_t, Grest_t, cos_sin_q.expMat, cos_sin_q.id);
    end


end