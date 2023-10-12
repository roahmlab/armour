function [pZ] = remove_dependence(pZ_in, max_id)
% remove dependence (i.e. make independent) all generators that depend on
% an id greater than the max_id
id_idx = (pZ_in.id > max_id);
gen_idx = any(pZ_in.expMat(id_idx, :) ~= 0, 1);
pZ = polyZonotope_ROAHM(pZ_in.c, pZ_in.G(:, ~gen_idx), [pZ_in.Grest, pZ_in.G(:, gen_idx)], pZ_in.expMat(~id_idx, ~gen_idx), pZ_in.id(~id_idx));
end

