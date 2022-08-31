function [ExpSO3_] = expSO3(x)
% Exponential map on SO(3)
%
% :param input: element of the lie algebra so(3), represented as a vector with 3 components.       
%
% :returns: element of the group SO(3), i.e. 3x3 rotation matrix
        
    a  = norm(x,2);
    
    if a == 0
        ExpSO3_ = eye(3);
    elseif a > 1e-20
        alpha = sin(a)/a;
        beta = (1-cos(a))/a^2;
        ExpSO3_ = eye(3) + alpha*hat(x) + beta*hat(x)^2;
    else 
        ExpSO3_ = zeros(3);
        powi = eye(3);
        for i=0:15
              ExpSO3_ = ExpSO3_...
                      +(1/factorial(i))*powi;
              powi = powi * hat(x);
        end     
    end

end
