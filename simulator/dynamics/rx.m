function R = rx(th)

    c = cos(th);
    s = sin(th);

%     R = [1  0  0  0;
%          0  c -s  0;
%          0  s  c  0; 
%          0  0  0  1];

    R = [1  0  0;
         0  c -s;
         0  s  c];   
end


