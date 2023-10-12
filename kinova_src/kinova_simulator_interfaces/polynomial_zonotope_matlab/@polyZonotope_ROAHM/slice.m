function [y] = slice(pz, x)
% 20220411
% assuming Grest is empty
% x is a vector corresponding to ids 1, 2, ..., max_id
% returns a double by evaluating generator coefficients and adding to center

x = x(pz.id);
if isempty(pz.G)
    y = pz.c;
else
    monom = prod(x.^pz.expMat, 1);
    terms = pz.G.*monom;
    y = pz.c + sum(terms, 2);
end

end

