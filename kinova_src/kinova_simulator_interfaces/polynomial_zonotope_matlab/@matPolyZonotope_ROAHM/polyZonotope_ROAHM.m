function [pz] = polyZonotope_ROAHM(mpz)
% converts matPolyZonotope_ROAHM to polyZonotope_ROAHM by stacking columns
% of matrices into vectors

c = mpz.C(:);
g = [];
grest = [];

for i = 1:size(mpz.G, 3)
   g_tmp = mpz.G(:, :, i);
   g = [g, g_tmp(:)];
end

for i = 1:size(mpz.Grest, 3)
   grest_tmp = mpz.Grest(:, :, i);
   grest = [grest, grest_tmp(:)];
end

pz = polyZonotope_ROAHM(c, g, grest, mpz.expMat, mpz.id);

end

