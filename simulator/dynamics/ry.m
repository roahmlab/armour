function R = ry(th)

    c = cos(th);
    s = sin(th);

    R = [c  0  s;
         0  1  0;
        -s  0  c];         
end


