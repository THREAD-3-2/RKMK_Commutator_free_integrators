function u = deriv2(v,vdot,v2dot)
% Support function used to define the control functions
%
% :param v: vector 
% :param vdot: time derivative of vector v  
% :param v2dot: second time derivative of vector v  
%
% :returns: time derivative of the output of deriv1

% recall deriv1:   u = (vdot*norma-v*normadot)/norma^2;

    norma = vecnorm(v);
    norma1d = (vdot'*v)/norma;
    norma2d = (norma*((vdot'*vdot)+(v'*v2dot))-norma1d*(v'*vdot))/norma^2;
    
    u = 1/norma^3*(-2*vdot*norma*norma1d + norma^2*v2dot - norma*v*norma2d +...
        2*v*norma1d^2);
    
%   u = (v2dot*norma^3-v*norma2d*norma^2-2*vdot*norma1d*norma^2+2*v*norma1d^2*norma)/norma^4;

end
