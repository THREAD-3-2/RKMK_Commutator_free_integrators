function u = deriv3(v,vdot,v2dot,v3dot)
% Support function used to define the control functions
%
% :param v: vector 
% :param vdot: time derivative of vector v  
% :param v2dot: second time derivative of vector v  
% :param v3dot: third time derivative of vector v  
%
% :returns: time derivative of the output of deriv2

    norma = vecnorm(v);
    norma1d = (vdot'*v)/norma;
    norma2d = (norma*((vdot'*vdot)+(v'*v2dot))-norma1d*(v'*vdot))/norma^2;
    helpnum2 = norma2d*norma^2;
    norma3d = ((norma1d*((v2dot'*v)+(vdot'*vdot))+norma*...
              ((v3dot'*v)+3*(v2dot'*vdot))-norma2d*(vdot'*v) - norma1d*...
              ((v2dot'*v)+(vdot'*vdot)))*norma^2 - helpnum2*2*norma*norma1d)/norma^4;     
    
    u = 1/norma^4*(-3*vdot*norma^2*norma2d+6*vdot*norma*norma1d^2+...
        norma^3*v3dot-norma^2*v*norma3d - 3*norma^2*norma1d*v2dot+...
        6*norma*v*norma1d*norma2d-6*v*norma1d^3);

    % recall deriv3u = (v2dot*norma^3-v*norma2d*norma^2-2*vdot*norma1d*norma^2+2*v*norma1d^2*norma)/norma^4;

%     u = (v3dot*norma^3 + v2dot*norma1d*norma^2 - 3*vdot*norma2d*norma^2 -...
%         v*norma3d*norma^2 + 2*v*norma2d*norma1d*norma - 2*vdot*norma1d^2*norma +...
%         2*v*norma1d^3)/norma^8;

end
