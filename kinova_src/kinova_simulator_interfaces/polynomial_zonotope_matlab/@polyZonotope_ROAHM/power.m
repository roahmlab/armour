function [pZout] = power(pZ, n)
% assuming pZ is 1-D, repeatedly multiply by itself n times

if (size(pZ.c, 1) ~= 1)
    error('pZ must be 1-dimensional');
end

if (n == 0)
    pZout = 1;
    return;
end

pZout = pZ;

if (n == 1)
    return;
end

for i = 2:n
    pZout = pZout.*pZ;
end

end

