function u=derivhatq2w(q,w,qdot,wdot)
     
    % u = (q'*w)*q - (q'*q)*w
    
    u = (qdot'*w + q'*wdot)*q + (q'*w)*qdot - 2*(q'*qdot)*w - (q'*q)*wdot;

end