function [B,dB,ddB] = Bezier_kernel_deg5(t,id,order)
% give values of Bernstein polynomials at s
% Input:
    % s should be a 1xn array, s is the phase variable of a Bezier curve, all elements should be inside [0,1]
    % id should be 1,2,3...,6. The function will then return the id th Bernstein polynomial
    % order should be 0,1,2. This indicates whether B, dB or ddB is the first output argument. This is prioritized than the number of output arguments
% Output:
    % B: Berntein polynomials at phase variable s, should be [6xn] array
    % B: 1st order derivative of Berntein polynomials at phase variable s, should be [6xn] array
    % B: 2nd-order derivative Berntein polynomials at phase variable s, should be [6xn] array

if nargin < 2
    order = [];
end

if isa(t, 'double') || isa(t, 'sym')
    t2 = t.*2.0;
    t3 = t.^2;
    t4 = t.^3;
    t6 = t-1.0;
    t8 = t.*1.2e+2;
    t15 = t.*6.0e+2;
    t16 = t.*1.2e+3;
    t5 = t3.^2;
    t9 = t2-2.0;
    t10 = t4.*2.0e+1;
    t11 = t3.*6.0e+1;
    t12 = t6.^2;
    t13 = t6.^3;
    t17 = t3.*1.8e+2;
    t7 = t5.*5.0;
    t14 = t12.^2;
    t19 = t13.*2.0e+1;
    t20 = t12.*6.0e+1;
    t24 = t12.*1.8e+2;
    t25 = t.*t9.*1.8e+2;
    t26 = t3.*t12.*3.0e+1;
    t18 = t14.*5.0;
    B = [-t6.^5;t.*t18;t3.*t13.*-1.0e+1;t4.*t12.*1.0e+1;t5.*t6.*-5.0;t.^5]';

    if nargin > 1
        B = B(id);
    end

    if nargout > 1 || (nargin > 2 && order == 1)
        t21 = t.*t19;
        dB = [-t18;t18+t21;-t26-t.*t13.*2.0e+1;t26+t4.*t9.*1.0e+1;-t7-t4.*t6.*2.0e+1;t7]';

        if nargin > 1
            dB = dB(id);
        end

        if order == 1
            B = dB;
        end
    end
    if nargout > 2 || (nargin > 2 && order == 2)
        t22 = t.*t20;
        t23 = -t19;
        ddB = [t23;t13.*4.0e+1+t22;t23-t.*t12.*1.2e+2-t3.*t9.*3.0e+1;t10+t22+t9.*t11;t4.*-4.0e+1-t3.*t6.*6.0e+1;t10]';

        if nargin > 1
            ddB = ddB(id);
        end

        if order == 2
            B = ddB;
        end
    end
elseif isa(t, 'interval')
    bezierDegree = 5;

    B = interval(zeros(1,bezierDegree+1));

    sinf = t.inf;
    ssup = t.sup;
   
    for i = 1:(bezierDegree+1)
        maxima = (i - 1) / bezierDegree; % this is the maxima of i th Bernstein polynomial
        
        if sinf >= maxima
            sup = Bezier_kernel_deg5(sinf,i);
            inf = Bezier_kernel_deg5(ssup,i);
        elseif ssup <= maxima
            sup = Bezier_kernel_deg5(ssup,i);
            inf = Bezier_kernel_deg5(sinf,i);
        else
            sup = Bezier_kernel_deg5(maxima,i);
            inf = min(Bezier_kernel_deg5(sinf,i), Bezier_kernel_deg10(ssup,i));
        end
    
        B(i) = interval(inf,sup);
    end

    if nargout > 1
        dB = interval(zeros(1,bezierDegree+1));

        dBezierKernel_maxima_minima = [0, 0;
            0, 0.4;
            2/5 - 6^(1/2)/10, 6^(1/2)/10 + 2/5;
            3/5 - 6^(1/2)/10, 6^(1/2)/10 + 3/5;
            0.6, 1;
            1, 1];
    
        for i = 1:(bezierDegree+1)
            maxima = dBezierKernel_maxima_minima(i,1); 
		    minima = dBezierKernel_maxima_minima(i,2);
    
		    if ssup <= maxima
			    sup = Bezier_kernel_deg5(ssup,i,1);
			    inf = Bezier_kernel_deg5(sinf,i,1);
		    elseif sinf < maxima
			    sup = Bezier_kernel_deg5(maxima,i,1);
			    inf = min(Bezier_kernel_deg5(sinf,i,1), Bezier_kernel_deg5(ssup,i,1));
		    elseif ssup <= minima
			    sup = Bezier_kernel_deg5(sinf,i,1);
			    inf = Bezier_kernel_deg5(ssup,i,1);
		    elseif sinf < minima
			    sup = max(Bezier_kernel_deg5(sinf,i,1), Bezier_kernel_deg5(ssup,i,1));
			    inf = Bezier_kernel_deg5(minima,i,1);
            else
			    sup = Bezier_kernel_deg5(ssup,i,1);
			    inf = Bezier_kernel_deg5(sinf,i,1);
            end
    
		    dB(i) = interval(inf, sup);
        end
    end
else
    error('unrecognized data type')
end