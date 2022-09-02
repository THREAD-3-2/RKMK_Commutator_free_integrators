function u = deriv1(v,vdot)
% Support function used to define the control functions
%
% :param v: vector  
% :param vdot: time derivative of vector v
%
% :returns: time derivative of the vector v divided by the norm of v

    norma = vecnorm(v);
    normadot = (vdot'*v)/norma;
    u = (vdot*norma-v*normadot)/norma^2;

end
