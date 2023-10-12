function [alpha] = bernstein_to_poly(beta, n)
% converts bernstein polynomial coefficients to
% monomial coefficients (power coefficients?)
    
    for i = 0:n
        alpha{i+1} = 0;
        for j = 0:i
            alpha{i+1} = alpha{i+1} + (-1)^(i-j)*nchoosek(n, i)*nchoosek(i, j)*beta{j+1};
        end
    end
end

