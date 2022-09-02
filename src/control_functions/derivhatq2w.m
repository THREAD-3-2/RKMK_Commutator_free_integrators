function u=derivhatq2w(q,w,qdot,wdot)
% Support function used to define the control functions

    % u = (q'*w)*q - (q'*q)*w
    
    u = (qdot'*w + q'*wdot)*q + (q'*w)*qdot - 2*(q'*qdot)*w - (q'*q)*wdot;

end
