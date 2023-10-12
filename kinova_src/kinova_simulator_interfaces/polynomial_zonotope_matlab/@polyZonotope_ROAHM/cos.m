function pZ = cos(pZ_inp, order)
% represent cos of a polynomial zonotope as an overapproximative polynomial zonotope.
% author: Bohao Zhang
% updated by: Patrick Holmes for lagrange remainder
% based on code from CORA_2021/contSet/@taylm/cos

if nargin == 1
    order = 6;
end

pZ_c = cos(pZ_inp.c);

pZ = pZ_c;

cs_cf = cos(pZ_inp.c);
sn_cf = sin(pZ_inp.c);

pZ_Grest = [];

factor = 1;
T_factor = 1;
pZ_neighbor = pZ_inp - pZ_inp.c;

for i = 1:order
    factor = factor * i;
    T_factor = T_factor .* pZ_neighbor;
    if mod(i,2) == 0
        pZ = pZ + ((sgn_cs(i) * cs_cf / factor) .* T_factor);
    else
        pZ = pZ + ((sgn_cs(i) * sn_cf / factor) .* T_factor);
    end
end

% add lagrange remainder interval to Grest
rem = interval(pZ_neighbor);
remPow = interval(T_factor.*pZ_neighbor);

if mod(order + 1, 2) == 1
    J0 = sin(pZ_inp.c + interval(0, 1).*rem);
else
    J0 = cos(pZ_inp.c + interval(0, 1).*rem);
end

if mod(order,4) == 0 || mod(order,4) == 1
    J = -J0;
else
    J = J0;
end

remainder = 1/(factor*(order+1))*remPow.*J;

pZ.c = pZ.c + center(remainder);
pZ.Grest = [pZ.Grest, rad(remainder)];
pZ.Grest = sum(abs(pZ.Grest),2);

% % add lagrange remainder interval to Grest
% if mod(order + 1, 2) == 1
%     J0 = @(x) sin(x);
% else
%     J0 = @(x) cos(x);
% end

% % sample f^(n+1)(x*)(x-x0)^n
% sampleNum = 200;
% sampleIntermediateNum = 5;
% sampleData = zeros(sampleNum,sampleIntermediateNum);
% for i = 1:sampleNum
%     sx = randPoint(pZ_inp);
%     intermediatePos = [0,unifrnd(0,1,1,sampleIntermediateNum-2),1];
%     for j = 1:sampleIntermediateNum
%         lambda = intermediatePos(j) * pZ_inp.c + (1 - intermediatePos(j)) * sx;
%         sampleData(i,j) = J0(lambda) * (sx - pZ_inp.c) ^ (order + 1);
%     end
% end
% remPowJ = interval(min(min(sampleData)),max(max(sampleData)));

% if mod(order,4) == 0 || mod(order,4) == 1
%     remPowJ = -remPowJ;
% end

% remainder = 1 / factorial(order + 1) * remPowJ;

% pZ.c = pZ.c + center(remainder);
% pZ.Grest = [pZ.Grest, rad(remainder)];
% pZ.Grest = sum(abs(pZ.Grest),2);

% pZ = reduce(pZ,'girard',reduceParam);
% pZ.Grest = sum(abs(pZ.Grest),2);

end

function res = sgn_cs(n)
    if mod(n,4) == 0 || mod(n,4) == 3
        res = 1;
    else
        res = -1;
    end
end