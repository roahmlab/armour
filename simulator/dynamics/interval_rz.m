function R = interval_rz(th)

    R_inf = rz(infimum(th));
    R_sup = rz(supremum(th));
    
    [m, n] = size(R_inf);
    
    R = interval(zeros(m,n), zeros(m,n));
   
    for i = 1:m
        for j =1:n
            v1 = R_inf(i,j);
            v2 = R_sup(i,j);
            if v2 > v1
                R(i,j) = interval(v1,v2);
            else 
                R(i,j) = interval(v2, v1);
            end
        end
    end
    
end

