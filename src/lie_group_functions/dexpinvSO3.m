function A = dexpinvSO3(v,input)
% Inverse of the derivative of the exponential map on SO(3)
%
% :param v: element in the lie algebra so(3), represented as a 3x1 vector   
% :param input: element in the lie algebra so(3), represented as a 3x1 vector   
%
% :returns: inverse of dexp_v(input) as 3x1 vector

    B = hat(v);
    alpha = norm(v,2);
    tol = 1e-20;
    
    if alpha>tol
        func = ( 1-alpha * cot(alpha/2)/2 )/(alpha^2);
        dexpinvB = eye(3) - 0.5 * B + func*B*B;
    else
        funclow = 1/12 + alpha^2/720 + alpha^4/(30240);
        dexpinvB = eye(3) - 0.5 * B + funclow*B*B;
    end
        
    A = dexpinvB * input;
    
end
