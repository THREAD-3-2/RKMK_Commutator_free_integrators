function u = deriv1(v,vdot)

    norma = vecnorm(v);
    normadot = (vdot'*v)/norma;
    u = (vdot*norma-v*normadot)/norma^2;

end