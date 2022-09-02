function uperp = getUperp(m1,m2,my,L1,L2,B,trajectory)
% Support function used to define the control functions and the vector field

% :returns: 6x1 vector, with uperp1=uperp(1:3) and uperp2=uperp(4:6)

% equation numbers are a reference to Ref[29] in the paper related to the code

    
    ky = 12;
    kydot = 5;
    kq = 16;
    kw = 5;
    
    g = 9.81;
    e3 = [0;0;1];

    yd = trajectory(:,1);
    yd_dot = trajectory(:,2);
    yd_ddot = trajectory(:,3);
    yd_dddot = trajectory(:,4);
    yd_ddddot = trajectory(:,5);

    r1d = [sin(deg2rad(30));0;cos(deg2rad(30))];
    r2d = [-sin(deg2rad(30));0;cos(deg2rad(30))];

    sd = -[0;1;0];


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
    
    q1dot = hat(w1)*q1;       % da (1)
    q2dot = hat(w2)*q2;       % da (1)

    %   eyddot = Mq(B(:,:,i))\h2(B(:,:,i))-yd_ddot(:,i);  %ho deciso di usare invece eq(13)
    yddot = (mu1+mu2)/my + g*e3;     %eq (13)
    eyddot = yddot-yd_ddot;

    %   recall  Fd = my*(-ky*ey - kydot*eydot + yd_ddot - g*[0;0;1]);  with ky, kydot const
    Fd_dot = my*(-ky*eydot - kydot*eyddot + yd_dddot);

    mu1dot = costmu1*((Fd_dot'*q1 + Fd'*q1dot)*q1+(Fd'*q1)*q1dot);     % derivata (19)
    mu2dot = costmu2*((Fd_dot'*q2 + Fd'*q2dot)*q2+(Fd'*q2)*q2dot);     % derivata (19)

    ydddot = (mu1dot+mu2dot)/my;         %eq (13)
    eydddot = ydddot - yd_dddot;      % derivata terza (15)
    Fd_ddot = my*(-ky*eyddot - kydot*eydddot + yd_ddddot);

    Q = [-hat(Fd)^2*sd/norm(hat(Fd)^2*sd),...
         -hat(Fd)*sd/norm(hat(Fd)*sd),...
         -Fd/norm(Fd)];       % da (17)

    %   calcoliamo Qdot
    vec1 = -hat(Fd)^2*sd;
    vec1dot = -derivhatq2w(Fd,sd,Fd_dot,zeros(3,1));
    vec2 = -hat(Fd)*sd;
    vec2dot = -hat(Fd_dot)*sd;
    vec3 = -Fd;
    vec3dot = -Fd_dot;

    Qdot = [deriv1(vec1,vec1dot),...
            deriv1(vec2,vec2dot),...
            deriv1(vec3,vec3dot)];

    %   calcoliamo Qddot
    %     vec1ddot = ((Fd*Fd_ddot')+(Fd_ddot*Fd')+2*(Fd_dot*Fd_dot'))*sd;
    vec1ddot = -deriv2hatq2w(Fd,sd,Fd_dot,zeros(3,1),Fd_ddot,zeros(3,1));
    vec2ddot = -hat(Fd_ddot)*sd;
    vec3ddot = -Fd_ddot;

    Qddot = [deriv2(vec1,vec1dot,vec1ddot),...
             deriv2(vec2,vec2dot,vec2ddot),...
             deriv2(vec3,vec3dot,vec3ddot)];


    q1d = Q*r1d;           % (18)
    q2d = Q*r2d;           % (18) 
    q1d_dot = Qdot*r1d;
    q2d_dot = Qdot*r2d;
    w1d = cross(q1d,q1d_dot);     % definito dopo (23)
    w2d = cross(q2d,q2d_dot);     %  "       "      " 
    eq1 = cross(q1d, q1);    % (23)
    eq2 = cross(q2d, q2);    % (23)
    ew1 = w1 + (hat(q1))^2*w1d;     % (23)
    ew2 = w2 + (hat(q2))^2*w2d;     % (23)

    q1d_ddot = Qddot*r1d;
    q2d_ddot = Qddot*r2d;
    w1d_dot = cross(q1d,q1d_ddot);
    w2d_dot = cross(q2d,q2d_ddot);

%   we can write now uperp1 and uperp2, u1 and u2, h5 and h6
    uperp1 = m1*L1*hat(q1)*(-kq*eq1 - kw*ew1 - ...
             (q1'*w1d)*q1dot - hat(q1)^2*w1d_dot) -...
             (m1/my)*hat(q1)^2*mu2;
         
    uperp2 = m2*L2*hat(q2)*(-kq*eq2 - kw*ew2 - ...
             (q2'*w2d)*q2dot - hat(q2)^2*w2d_dot) -...
             (m2/my)*hat(q2)^2*mu1;
     
     uperp = [uperp1;uperp2];
     
         
end
