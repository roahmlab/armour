function R = rz(th)

    c = cos(th);
    s = sin(th);

%     R = [c -s 0  0;
%          s  c 0  0; 
%          0  0 1  0;
%          0  0 0  1];

    R = [c -s 0;
         s  c 0; 
         0  0 1];
end