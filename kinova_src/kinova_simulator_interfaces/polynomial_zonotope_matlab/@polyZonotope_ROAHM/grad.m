function [grad] = grad(pz,max_id)
% 20220412
% returns the gradient of the pz with respect to id 1, 2, ... max_id
% outputs a max_id x 1 cell of poly zonotopes, which may be sliced
% to get the gradient

grad = cell(max_id, 1);

for i = 1:max_id
    c_tmp = zeros(size(pz.c));
    idx = find(pz.id == i, 1);
    if isempty(idx)
       grad{i, 1} = polyZonotope_ROAHM(zeros(size(pz.c)));
       continue;
    end
    G_tmp = pz.expMat(idx, :).*pz.G;
    E_tmp = pz.expMat;
    E_tmp(idx, :) = E_tmp(idx, :) - 1;
    E_tmp(E_tmp < 0) = 0;

    empty_idx = all(E_tmp == 0, 1);
    c_tmp = c_tmp + sum(G_tmp(:, empty_idx), 2);
    G_tmp(:, empty_idx) = [];
    E_tmp(:, empty_idx) = [];


    grad{i, 1} = polyZonotope_ROAHM(c_tmp, G_tmp, [], E_tmp, pz.id);
end


end

