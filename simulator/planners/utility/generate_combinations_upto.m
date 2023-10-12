function [combs] = generate_combinations_upto(N)
% generate a bunch of combinations, store in a cell

combs = cell(N, 1);
combs{1} = [1];

for i = 2:N
    combs{i} = combinator(i, 2, 'c');
end

end

