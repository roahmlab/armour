function T = create_interval(T1, T2)

[m,n] = size(T1);

Tinf = zeros(m,n);
Tsup = zeros(m,n);
for i = 1:m
    for j = 1:n
        a = min(T1(i,j), T2(i,j));
        b = max(T1(i,j), T2(i,j));
        
        Tinf(i,j) = a;
        Tsup(i,j) = b;
    end
end

T = interval(Tinf, Tsup);

end

