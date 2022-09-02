function u=deriv3hatq2w(q,w,qdot,wdot,q2dot,w2dot,q3dot,w3dot)
% Support function used to define the control functions

%   u = (q'*w)*q - (q'*q)*w;

%   u = (qdot'*w + q'*wdot)*q + (q'*w)*qdot - 2*(q'*qdot)*w - (q'*q)*wdot;

%   u = (q2dot'*w + 2*qdot'*wdot + q'*w2dot)*q + 2*(qdot'*w + q'*wdot)*qdot +...
%       (q'*w)*q2dot - 2*(qdot'*qdot + q'*q2dot)*w - 4*(q'*qdot)*wdot -...
%       (q'*q)*w2dot;
    
    u = (q3dot'*w + 3*q2dot'*wdot + 3*qdot'*w2dot + q'*w3dot)*q +...
        3*(q2dot'*w + 2*qdot'*wdot + q'*w2dot)*qdot + 3*(qdot'*w + q'*wdot)*q2dot +...
        (q'*w)*q3dot -2*(3*qdot'*q2dot + q'*q3dot)*w - 6*(qdot'*qdot + q'*q2dot)*wdot -...
        6*(q'*qdot)*w2dot - (q'*q)*w3dot;
              
    
end
