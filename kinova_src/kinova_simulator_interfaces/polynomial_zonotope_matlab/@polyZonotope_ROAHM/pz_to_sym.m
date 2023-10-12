function [out] = pz_to_sym(pz, x)
% converts polyZonotope_ROAHM to a symbolic expression
% assumes that Grest is empty, and that all generators will be fully sliced
% (i.e. all ids are sliced at a point)
% requires that the symbolic vector x is passed in

if ~isempty(pz.Grest)
   error('Assumes Grest is empty (only dependent generators)'); 
end

x = x(pz.id);
monom = prod(x.^pz.expMat, 1);
terms = pz.G.*monom;
out = pz.c + sum(terms, 2);

end

