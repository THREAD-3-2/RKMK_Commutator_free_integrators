function upar = getUpar(m1,m2,my,L1,L2,B,trajectory)
% Support function used to define the control functions and the vector field

% :returns: 6x1 vector, with upar(1:3)=upar1 and upar(4:6)=upar2

% equation numbers are a reference to Ref[29] in the paper related to the code

    ky = 12;
    kydot = 5;
    
    g = 9.81;
    e3 = [0;0;1];



    yd = trajectory(:,1);
    yd_dot = trajectory(:,2);
    yd_ddot = trajectory(:,3);

    r1d = [sin(deg2rad(30));0;cos(deg2rad(30))];
    r2d = [-sin(deg2rad(30));0;cos(deg2rad(30))];

    q1 = B(:,11);    
    q2 = B(:,13);    
    w1 = B(:,12);    
    w2 = B(:,14);    

    y = B(:,1);
    ydot = B(:,2);

%   definition upar

    ey = y - yd;     % eq (15)
    eydot = ydot - yd_dot;
    Fd = my*(-ky*ey - kydot*eydot + yd_ddot - g*[0;0;1]);    % eq (16)
    costmu1 = 1/((r1d'*e3)*norm(r1d+r2d));
    mu1 = costmu1*(Fd'*q1)*q1;   % eq (19)
    costmu2 = 1/((r2d'*e3)*norm(r1d+r2d));
    mu2 = costmu2*(Fd'*q2)*q2;   % eq (19)

    upar1 = m1*L1*(norm(w1))^2*q1 + mu1 +...        
            (m1/my)*(q1*q1')*(mu1+mu2);                % eq (12)
    upar2 = m2*L2*(norm(w2))^2*q2 + mu2 +...
            (m2/my)*(q2*q2')*(mu1+mu2);
        
    upar = [upar1;upar2];
end
