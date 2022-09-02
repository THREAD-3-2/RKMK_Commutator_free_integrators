function rc = controls(m1,m2,my,L1,L2,J1,J2,B,trajectory)
% Main function used to define the control functions

% equation numbers are a reference to Ref[29] in the paper related to the code


    ky = 12;
    kydot = 5;
    kq = 16;
    kw = 5;
    epsi = 0.9;
    kR = 9;
    kW = 3;
    
    g = 9.81;
    e3 = [0;0;1];
    
    yd = trajectory(:,1);
    yd_dot = trajectory(:,2);
    yd_ddot = trajectory(:,3);
    yd_dddot = trajectory(:,4);
    yd_ddddot = trajectory(:,5);
    yd_dddddot = trajectory(:,6);
    yd_ddddddot = trajectory(:,7);

    r1d = [sin(deg2rad(30));0;cos(deg2rad(30))];
    r2d = [-sin(deg2rad(30));0;cos(deg2rad(30))];

    sd = -[0;1;0];
    b11 = -[1;0;0];
    b12 = -[1;0;0];


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

    q1dot = hat(w1)*q1;       % da (1)
    q2dot = hat(w2)*q2;       % da (1)

    %   eyddot = Mq(B(:,:,i))\h2(B(:,:,i))-yd_ddot(:,i);  %instead of this we used eq(13)
    yddot = (mu1+mu2)/my + g*e3;     %eq (13)
    eyddot = yddot-yd_ddot;

    %   recall  Fd = my*(-ky*ey - kydot*eydot + yd_ddot - g*[0;0;1]);  with ky, kydot const
    Fd_dot = my*(-ky*eydot - kydot*eyddot + yd_dddot);

    mu1dot = costmu1*((Fd_dot'*q1 + Fd'*q1dot)*q1+(Fd'*q1)*q1dot);     % derivative (19)
    mu2dot = costmu2*((Fd_dot'*q2 + Fd'*q2dot)*q2+(Fd'*q2)*q2dot);     % derivative (19)

    ydddot = (mu1dot+mu2dot)/my;         %eq (13)
    eydddot = ydddot - yd_dddot;      % derivata terza (15)
    Fd_ddot = my*(-ky*eyddot - kydot*eydddot + yd_ddddot);

    Q = [-hat(Fd)^2*sd/norm(hat(Fd)^2*sd),...
         -hat(Fd)*sd/norm(hat(Fd)*sd),...
         -Fd/norm(Fd)];       % da (17)

    %   let us evaluate Qdot
    vec1 = -hat(Fd)^2*sd;
    vec1dot = -derivhatq2w(Fd,sd,Fd_dot,zeros(3,1));
    vec2 = -hat(Fd)*sd;
    vec2dot = -hat(Fd_dot)*sd;
    vec3 = -Fd;
    vec3dot = -Fd_dot;

    Qdot = [deriv1(vec1,vec1dot),...
            deriv1(vec2,vec2dot),...
            deriv1(vec3,vec3dot)];

    %   let us evaluate Qddot
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
    w1d = cross(q1d,q1d_dot);     % defined after (23)
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

    u1 = upar1 + uperp1;
    u2 = upar2 + uperp2;

    
    b31 = - u1/norm(u1);        % (27)
    b32 = - u2/norm(u2);        % (27)
    Rc1 = [-hat(b31)^2*b11/norm(hat(b31)^2*b11) , ...
            hat(b31)*b11/norm(hat(b31)*b11) , b31];   %(28)
    Rc2 = [-hat(b32)^2*b12/norm(hat(b32)^2*b12) , ...
            hat(b32)*b12/norm(hat(b32)*b12) , b32];   %(28)
   
    
%     w1dot= 1/(my*L1)*hat(sol(:,11,i))*mu2 - ...
%               1/(m1*L1)*hat(sol(:,11,i))*uperp1;    %eq (14)
%     w2dot = 1/(my*L2)*hat(sol(:,13,i))*mu1 - ...
%               1/(m2*L2)*hat(sol(:,13,i))*uperp2;
          
    w1dot = 1/L1*hat(q1)*(yddot-g*e3)-1/(m1*L1)*hat(q1)*uperp1;   % eq (6)
    w2dot = 1/L2*hat(q2)*(yddot-g*e3)-1/(m2*L2)*hat(q2)*uperp2;

          
    %    derivative of upar (12)
    upar1_dot = 2*m1*L1*(w1'*w1dot)*q1 + m1*L1*norm(w1)^2*q1dot+mu1dot +...
                (m1/my)*( ((q1dot*q1') + (q1*q1dot'))*(mu1+mu2) +...
                ((q1*q1')*(mu1dot+mu2dot)) );
    upar2_dot = 2*m2*L2*(w2'*w2dot)*q2 + m2*L2*norm(w2)^2*q2dot+mu2dot +...
                (m2/my)*(((q2dot*q2') + (q2*q2dot'))*(mu1+mu2) +...
                ((q2*q2')*(mu1dot+mu2dot)));
        
    % now we need the derivative of uperp 
    eq1dot = cross(q1d_dot,q1)+cross(q1d,q1dot);      % cf (23)
    eq2dot = cross(q2d_dot,q2)+cross(q2d,q2dot);
    
%     w1dotbis = -kq*eq1-kw*ew1-(q1*w1d')*q1dot-hat(q1)^2*w1d_dot;    % (24)
%     w2dotbis = -kq*eq2-kw*ew2-(q2*w2d')*q2dot-hat(q2)^2*w2d_dot;    % (24) 
    
%     ew1dot = w1dot + hat(sol(:,11,i))^2*w1dot - (sol(:,11,i)*q1dot'+q1dot*sol(:,11,i)')*w1d;
%     ew2dot = w2dot + hat(sol(:,13,i))^2*w2dot - (sol(:,13,i)*q2dot'+q2dot*sol(:,13,i)')*w2d;
       
    ew1dot = w1dot + derivhatq2w(q1,w1d,q1dot,w1d_dot);     % cf (23) 
    ew2dot = w2dot + derivhatq2w(q1,w2d,q2dot,w2d_dot);
    
    q1ddot = cross(w1dot,q1)+cross(w1,q1dot);      % cf (1)
    q2ddot = cross(w2dot,q2)+cross(w2,q2dot);
    
    %    Recall: mu1dot = costmu1*((Fd_dot'*sol(:,11,i) +...
    %                     Fd'*q1dot)*sol(:,11,i)+(Fd'*sol(:,11,i))*q1dot);     % derivata (19)
    mu1ddot = costmu1*((Fd_ddot'*q1 + 2*Fd_dot'*q1dot + Fd'*q1ddot)*q1 +...
              2*(Fd_dot'*q1+Fd'*q1dot)*q1dot + (Fd'*q1)*q1ddot);
    mu2ddot = costmu2*((Fd_ddot'*q2 + 2*Fd_dot'*q2dot + Fd'*q2ddot)*q2 +...
              2*(Fd_dot'*q2+Fd'*q2dot)*q2dot + (Fd'*q2)*q2ddot);
    
    % let us consider again the derivatives of Q
    
    y_ddddot = (mu1ddot+mu2ddot)/my;         %eq (13)
    eyddddot = y_ddddot - yd_ddddot;
    Fd_dddot = my*(-ky*eydddot - kydot*eyddddot + yd_dddddot);
    
%     vec1dddot = (Fd*Fd_dddot'+Fd_dddot*Fd'+3*Fd_ddot*Fd_dot'+...
%                  3*Fd_dot*Fd_ddot')*sd;
    vec1dddot = -deriv3hatq2w(Fd,sd,Fd_dot,zeros(3,1),Fd_ddot,zeros(3,1),Fd_dddot,zeros(3,1));    
    vec2dddot = -hat(Fd_dddot)*sd;
    vec3dddot = -Fd_dddot;
    
    Qdddot = [deriv3(vec1,vec1dot,vec1ddot,vec1dddot),...
              deriv3(vec2,vec2dot,vec2ddot,vec2dddot),...
              deriv3(vec3,vec3dot,vec3ddot,vec3dddot)];

    
    q1d_dddot = Qdddot*r1d;       % cf (18)
    q2d_dddot = Qdddot*r2d;
    
    w1d_ddot = cross(q1d_dot,q1d_ddot)+cross(q1d,q1d_dddot);     % cf after (23)
    w2d_ddot = cross(q2d_dot,q2d_ddot)+cross(q2d,q2d_dddot);
       
    % support variables for the derivatives of uperp1 and uperp2
    
    uperp1help1 = (-kq*eq1 - kw*ew1 - (q1'*w1d)*q1dot - hat(q1)^2*w1d_dot);
    uperp2help1 = (-kq*eq2 - kw*ew2 - (q2'*w2d)*q2dot - hat(q2)^2*w2d_dot);
    
    % %     uperp1help2 = hat(q1)^2*w1d_ddot-(q1*q1dot'+...
    % %                   q1dot*q1')*w1d_dot;
    % %     uperp2help2 = hat(q2)^2*w2d_ddot-(q2*q2dot'+...
    % %                   q2dot*q2')*w2d_dot;
    uperp1help2 = derivhatq2w(q1,w1d_dot,q1dot,w1d_ddot);
    uperp2help2 = derivhatq2w(q2,w2d_dot,q2dot,w2d_ddot);
    
    % %     dhq1mu = -m1/my*(hat(sol(:,11,i))^2*mu2dot -...
    % %               ((sol(:,11,i)*q1dot')+(q1dot*(sol(:,11,i))'))*mu2);  % d/dt (-mi/my)*hat(qi)^2*muj
    % %     dhq2mu = -m2/my*(hat(sol(:,13,i))^2*mu1dot -...
    % %               ((sol(:,13,i)*q2dot')+(q2dot*(sol(:,13,i))'))*mu1);
    uperp1help3 = -(m1/my)*derivhatq2w(q1,mu2,q1dot,mu2dot);
    uperp2help3 = -(m2/my)*derivhatq2w(q2,mu1,q2dot,mu1dot);
              
%   derivative of uperp1 and uperp2

    uperp1_dot = m1*L1*hat(q1dot)*uperp1help1 + m1*L1*hat(q1)*...
                 (-kq*eq1dot-kw*ew1dot - (q1dot'*w1d+q1'*w1d_dot)*q1dot -...
                 (q1'*w1d)*q1ddot- uperp1help2) + uperp1help3;
    uperp2_dot = m2*L2*hat(q2dot)*uperp2help1 + m2*L2*hat(q2)*...
                 (-kq*eq2dot-kw*ew2dot - (q2dot'*w2d+q2'*w2d_dot)*q2dot -...
                 (q2'*w2d)*q2ddot- uperp2help2) + uperp2help3;         

%   derivative of u1 and u2             
    u1dot = upar1_dot + uperp1_dot;
    u2dot = upar2_dot + uperp2_dot;
    
    % we have to write Rc1dot e Rc2dot derivatives of (28) with b3i given by (27)
%     b31_dot = -(u1dot/vecnorm(u1) - u1*((u1dot'*u1)/vecnorm(u1)^3));
%     b32_dot = -(u2dot/vecnorm(u2) - u2*((u2dot'*u2)/vecnorm(u2)^3));
    b31_dot = -deriv1(u1,u1dot);
    b32_dot = -deriv1(u2,u2dot);

    vec1b1 = -hat(b31)^2*b11;
        %     vec1b1_dot = (b31*b31_dot' + b31_dot*b31')*b11;   % derivative of -hat(b31)^2*b11
    vec1b1_dot = -derivhatq2w(b31,b11,b31_dot,zeros(3,1));
    vec1b2 = -hat(b32)^2*b12;
        %     vec1b2_dot = (b32*b32_dot' + b32_dot*b32')*b11;   % derivative of -hat(b32)^2*b11
    vec1b2_dot = -derivhatq2w(b32,b12,b32_dot,zeros(3,1));
    vec2b1 = hat(b31)*b11;
    vec2b1_dot = hat(b31_dot)*b11;
    vec2b2 = hat(b32)*b12;
    vec2b2_dot = hat(b32_dot)*b12;
    vec3b1 = -u1;
    vec3b2 = -u2;
    vec3b1_dot = -u1dot;
    vec3b2_dot = -u2dot;

    Rc1dot = [deriv1(vec1b1,vec1b1_dot),...
              deriv1(vec2b1,vec2b1_dot),...
              deriv1(vec3b1,vec3b1_dot)];
    
    Rc2dot = [deriv1(vec1b2,vec1b2_dot),...
              deriv1(vec2b2,vec2b2_dot),...
              deriv1(vec3b2,vec3b2_dot)];
          
    % other pieces of code for M1 e M2 cf. (30)
    Wc1 = invhat(Rc1'*Rc1dot);
    Wc2 = invhat(Rc2'*Rc2dot);
    eR1 = 0.5*invhat(Rc1'*B(:,3:5)-(B(:,3:5))'*Rc1);
    eR2 = 0.5*invhat(Rc2'*B(:,7:9)-(B(:,7:9))'*Rc2);
    eW1 = B(:,6)-(B(:,3:5))'*Rc1*Wc1;
    eW2 = B(:,10)-(B(:,7:9))'*Rc2*Wc2;
    
    
    % finally, for M1 and M2, we need Wc1_dot and Wc2_dot
    
    
    w1ddot = 1/(my*L1)*(hat(q1dot)*mu2 + hat(q1)*mu2dot) +...
             1/(m1*L1)*(hat(q1dot)*uperp1 + hat(q1)*uperp1_dot);      % cf (14)
    w2ddot = 1/(my*L2)*(hat(q2dot)*mu1 + hat(q2)*mu1dot) +...
             1/(m1*L2)*(hat(q2dot)*uperp2 + hat(q2)*uperp2_dot);
    
    % Recall:    q1ddot = cross(w1dot,q1)+cross(sol(:,12,i),q1dot);      % cf (1)
    q1dddot = cross(w1ddot,q1) + 2*cross(w1dot,q1dot) + cross(B(:,12),q1ddot);
    q2dddot = cross(w2ddot,q2) + 2*cross(w2dot,q2dot) + cross(B(:,14),q2ddot);
    
    % Recall:   mu1ddot = costmu*((Fd_ddot'*q1+2*Fd_dot'*q1dot+...
    %                     Fd'*q1ddot)*q1 + 2*(Fd_dot'*q1+Fd'*q1dot)*q1dot + ...
    %                    (Fd'*q1)*q1ddot);
    mu1dddhelp1 = (Fd_ddot'*q1 + 2*Fd_dot'*q1dot + Fd'*q1ddot);
    mu1dddhelp2 = (Fd_dot'*q1 + Fd'*q1dot);
    mu2dddhelp1 = (Fd_ddot'*q2 + 2*Fd_dot'*q2dot + Fd'*q2ddot);
    mu2dddhelp2 = (Fd_dot'*q2 + Fd'*q2dot);
    
    mu1dddot = costmu1*((Fd_dddot'*q1 + 3*Fd_ddot'*q1dot +...
               3*Fd_dot'*q1ddot + Fd'*q1dddot)*q1 +...
               3*mu1dddhelp1*q1dot + 3*mu1dddhelp2*q1ddot + (Fd'*q1)*q1dddot);
    mu2dddot = costmu2*((Fd_dddot'*q2 + 3*Fd_ddot'*q2dot +...
               3*Fd_dot'*q2ddot + Fd'*q2dddot)*q2 +...
               3*mu2dddhelp1*q2dot + 3*mu2dddhelp2*q2ddot + (Fd'*q2)*q2dddot);
    
    y_dddddot = (mu1dddot+mu2dddot)/my;         %da eq (13)
    eydddddot = y_dddddot - yd_dddddot;
    Fd_ddddot = my*(-ky*eyddddot - kydot*eydddddot + yd_ddddddot);
    
% %     vec1ddddot = (Fd*Fd_ddddot'+Fd_ddddot*Fd'+4*Fd_dddot*Fd_dot'+...
% %                   4*Fd_dot*Fd_dddot'+6*(Fd_ddot*Fd_ddot'))*sd;
    vec1ddddot = -deriv4hatq2w(Fd,sd,Fd_dot,zeros(3,1),Fd_ddot,zeros(3,1),Fd_dddot,zeros(3,1),Fd_ddddot,zeros(3,1));
    vec2ddddot = -hat(Fd_ddddot)*sd;
    vec3ddddot = -Fd_ddddot;
    
    Qddddot = [deriv4(vec1,vec1dot,vec1ddot,vec1dddot,vec1ddddot),...
               deriv4(vec2,vec2dot,vec2ddot,vec2dddot,vec1ddddot),...
               deriv4(vec3,vec3dot,vec3ddot,vec3dddot,vec1ddddot)];      
    

    upar1_ddot = m1*L1*(2*(w1ddot'*w1 + vecnorm(w1dot)^2)*q1 +...
                 3*w1dot'*w1*q1dot + vecnorm(w1dot)^2*q1ddot) + ...
                 mu1ddot + m1/my*((q1ddot*q1' + 2*(q1dot*q1dot')+...
                 q1*q1ddot')*(mu1+mu2) + 2*(q1dot*q1' +...
                 q1*q1dot')*(mu1dot+mu2dot) +...
                 q1*q1'*(mu1ddot+mu2ddot));
    upar2_ddot = m2*L2*(2*(w2ddot'*w2 + vecnorm(w2dot)^2)*q2 +...
                 3*w2dot'*w2*q2dot + vecnorm(w2dot)^2*q2ddot) + ...
                 mu2ddot + m2/my*((q2ddot*q2' + 2*(q2dot*q2dot') +...
                 q2*q2ddot')*(mu1+mu2) + 2*(q2dot*q2' +...
                 q2*q2dot')*(mu1dot+mu2dot) +...
                 q2*q2'*(mu1ddot+mu2ddot));
      
    q1d_dddot = Qdddot*r1d;
    q2d_dddot = Qdddot*r2d;
    q1d_ddddot = Qddddot*r1d;
    q2d_ddddot = Qddddot*r2d;
    w1d_dddot = 2*cross(q1d_dot,q1d_dddot) + cross(q1d,q1d_ddddot);
    w2d_dddot = 2*cross(q2d_dot,q2d_dddot) + cross(q2d,q2d_ddddot);
    
    eq1ddot = cross(q1d_ddot,q1) + 2*cross(q1d_dot,q1dot) +...
              cross(q1d,q1ddot);
    eq2ddot = cross(q2d_ddot,q2) + 2*cross(q2d_dot,q2dot) +...
              cross(q2d,q2ddot);
    
    ew1ddot = w1ddot + deriv2hatq2w(q1,w1d,q1dot,w1d_dot,q1ddot,w1d_ddot);
    ew2ddot = w2ddot + deriv2hatq2w(q2,w2d,q2dot,w2d_dot,q2ddot,w2d_ddot);
% %     ew2ddot = w2ddot + hat(q2)^2*w2d_ddot-2*(q2*q2dot'+...
% %               q2dot*q2')*w2d_dot-(q2*q2ddot'+2*(q2dot*q2dot')+...
% %               q2ddot*q2')*w2d;
% %     termine4_1 = hat(q1)^2*w1d_dddot-2*(q1*q1dot'+...
% %                  q1dot*q1')*w2d_ddot-(q1*q1ddot'+2*(q1dot*q1dot')+...
% %                  q1ddot*q1')*w2d_dot;
    termine3_1 = (q1'*w1d)*q1dddot+2*(q1dot'*w1d+q1'*w1d_dot)*q1ddot+...
                 (q1ddot'*w1d+2*q1dot'*w1d_dot+q1'*w1d_ddot)*q1dot;
    termine4_1 = deriv2hatq2w(q1,w1d_dot,q1dot,w1d_ddot,q1ddot,w1d_dddot);
% %     termine4_2 = hat(q2)^2*w2d_dddot-2*(q2*q2dot'+...
% %                  q2dot*q2')*w2d_ddot-(q2*q2ddot'+2*(q2dot*q2dot')+...
% %                  q2ddot*q2')*w2d_dot;
% (q1dot'*w1d + sol'*w1dot)*q1dot + (sol'*w1d)*q1ddot
% 
    termine3_2 = (q2'*w2d)*q2dddot+2*(q2dot'*w2d+q2'*w2d_dot)*q2ddot+...
                 (q2ddot'*w2d+2*q2dot'*w2d_dot+q2'*w2d_ddot)*q2dot;
    termine4_2 = deriv2hatq2w(q2,w2d_dot,q2dot,w2d_ddot,q2ddot,w2d_dddot);
    
    help1a = (-kq*eq1dot-kw*ew1dot - (q1dot'*w1d+q1'*w1d_dot)*q1dot-...
             (q1'*w1d)*q1ddot- uperp1help2);
    help2a = (-kq*eq2dot-kw*ew2dot - (q2dot'*w2d+q2'*w2d_dot)*q2dot-...
             (q2'*w2d)*q2ddot- uperp2help2);
         
    help1b = -kq*eq1ddot-kw*ew1ddot-termine3_1-termine4_1;
    help2b = -kq*eq2ddot-kw*ew2ddot-termine3_2-termine4_2;
    
%     uperp1_ddot = m1*L1*hat(q1ddot)*uperp1help + 2*m1*L1*hat(q1dot)*help1a+...
%                   m1*L1*hat(q1)*help1b + m1/my*(hat(q1)^2*mu2ddot-...
%                   2*(q1*q1dot'+q1dot*q1')*mu2dot - ...
%                   (q1*q1ddot'+2*(q1dot*q1dot')+q1ddot*q1')*mu2);
%     uperp2_ddot = m2*L2*hat(q2ddot)*uperp2help + 2*m2*L2*hat(q2dot)*help2a+...
%                   m2*L2*hat(q2)*help2b + m2/my*(hat(q2)^2*mu1ddot-...
%                   2*(q2*q2dot'+q2dot*q2')*mu1dot - ...
%                   (q2*q2ddot'+2*(q2dot*q2dot')+q2ddot*q2')*mu1);       
%    
    uperp1_ddot = m1*L1*hat(q1ddot)*uperp1help1 + 2*m1*L1*hat(q1dot)*help1a+...
                  m1*L1*hat(q1)*help1b -...
                  m1/my*deriv2hatq2w(q1,mu2,q1dot,mu2dot,q1ddot,mu2ddot);
    uperp2_ddot = m2*L2*hat(q2ddot)*uperp2help1 + 2*m2*L2*hat(q2dot)*help2a+...
                  m2*L2*hat(q2)*help2b -...
                  m2/my*deriv2hatq2w(q2,mu1,q2dot,mu1dot,q2ddot,mu1ddot);


    u1ddot = upar1_ddot + uperp1_ddot;
    u2ddot = upar2_ddot + uperp2_ddot;
    
    b31_ddot = -deriv2(u1,u1dot,u1ddot);
    b32_ddot = -deriv2(u2,u2dot,u2ddot);

% %     vec1b1_ddot = (b31*b31_ddot' + b31_ddot*b31' + 2*(b31_dot*b31_dot'))*b11; 
% %     vec1b2_ddot = (b32*b32_ddot' + b32_ddot*b32' + 2*(b32_dot*b32_dot'))*b11;  
    vec1b1_ddot = -deriv2hatq2w(b31,b11,b31_dot,zeros(3,1),b31_ddot,zeros(3,1));
    vec1b2_ddot = -deriv2hatq2w(b31,b11,b31_dot,zeros(3,1),b32_ddot,zeros(3,1));
    vec2b1_ddot = hat(b31_ddot)*b11;
    vec2b2_ddot = hat(b32_ddot)*b11;
    vec3b1_ddot = -u1ddot;
    vec3b2_ddot = -u2ddot;

    Rc1ddot = [deriv2(vec1b1,vec1b1_dot,vec1b1_ddot),...
               deriv2(vec2b1,vec2b1_dot,vec2b1_ddot),...
               deriv2(vec3b1,vec3b1_dot,vec3b1_ddot)];
    Rc2ddot = [deriv2(vec1b2,vec1b2_dot,vec1b2_ddot),...
               deriv2(vec2b2,vec2b2_dot,vec2b2_ddot),...
               deriv2(vec3b2,vec3b2_dot,vec3b2_ddot)]; 
    
%     Wc1dot = invhat(Rc1dot'*Rc1dot+Rc1'*Rc1ddot);    
%     Wc2dot = invhat(Rc2dot'*Rc2dot+Rc2'*Rc2ddot);
    Wc1dot = invhat(Rc1'*(Rc1ddot-Rc1dot*hat(Wc1)));
    Wc2dot = invhat(Rc2'*(Rc2ddot-Rc2dot*hat(Wc2)));
    
%   finally M1 and M2, h3 and h4
    
    M1 = -(kR/epsi^2)*eR1 - (kW/epsi)*eW1 + hat(B(:,6))*J1*B(:,6) - ...
         J1*(hat(B(:,6))*(B(:,3:5))'*Rc1*Wc1 - (B(:,3:5))'*Rc1*Wc1dot);  % eq (30)
    M2 = -(kR/epsi^2)*eR2 - (kW/epsi)*eW2 + hat(B(:,10))*J2*B(:,10) - ...
         J2*(hat(B(:,10))*(B(:,7:9))'*Rc2*Wc2 - (B(:,7:9))'*Rc2*Wc2dot); % eq (30)
%     M1=M1/1000;
%     M2=M2/1000;
    
%     M1 = [0.1;0.02;0.01];
%     M2 = [0;0.22;0.07];
%     
%     
%     uperp1 = (eye(3)-q1*q1')*rand(3,1);
%     uperp2 = (eye(3)-q2*q2')*rand(3,1);

    rc = [q1d,q2d,Rc1,Rc2];

end
