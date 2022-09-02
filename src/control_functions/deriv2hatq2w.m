function u=deriv2hatq2w(q,w,qdot,wdot,q2dot,w2dot)
% Support function used to define the control functions

%   u = (q'*w)*q - (q'*q)*w;

%   u = (qdot'*w + q'*wdot)*q + (q'*w)*qdot - 2*(q'*qdot)*w - (q'*q)*wdot;

    u = (q2dot'*w + 2*qdot'*wdot + q'*w2dot)*q + 2*(qdot'*w + q'*wdot)*qdot +...
        (q'*w)*q2dot - 2*(qdot'*qdot + q'*q2dot)*w -4*(q'*qdot)*wdot -...
        (q'*q)*w2dot;

    
end
