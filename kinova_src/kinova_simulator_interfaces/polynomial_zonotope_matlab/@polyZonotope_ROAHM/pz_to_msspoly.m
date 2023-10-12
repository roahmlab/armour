function [out, x] = pz_to_msspoly(pz, x)
% converts polyZonotope_ROAHM to msspoly
    n_id = length(pz.id);
    n_gen = size(pz.G, 2);
    n_dim = size(pz.c, 1);
    
    %SHAN: added this. In levi_data_FRS_polyZonotope, we only input pz to
    %this function. Was throwing an error.
    if nargin == 1
        x = msspoly('x', n_id);
    end

    % reorder according to pz id order:
    x = x(pz.id);
    sz = n_gen + 1;
    dim = [n_dim, sz];
    sub = [];
    for i = 1:n_dim
        sub = [sub; [i*ones(sz, 1), (1:sz)']];
    end
    var = repmat([x.var'], [n_dim*sz 1]);
    pow = repmat([zeros(1, n_id); pz.expMat'], [n_dim 1]);
    coeff = [];
    for i = 1:n_dim
        if pz.c(i, 1) == 0
            % use a very small number
            coeff = [coeff; [realmin('double'); pz.G(i, :)']];
        elseif isempty(pz.G)
            coeff = [coeff; pz.c(i, 1)'];
        else
            coeff = [coeff; [pz.c(i, 1)'; pz.G(i, :)']];
        end
    end
    out = msspoly(dim, sub, var, pow, coeff);
    out = sum(out, 2);
end

