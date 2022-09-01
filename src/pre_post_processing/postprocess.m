% Postprocessing: update of support variables to plot errors and real time plot of the tracking problem

eynorm(i) = norm(sol(:,1)-yd(t),2);

q1plot(:,i+1) = sol(:,11);
q2plot(:,i+1) = sol(:,13);

y1plot(:,i+1) = sol(1,1);
y2plot(:,i+1) = sol(2,1);
y3plot(:,i+1) = sol(3,1);

mass(:,i+1) = sol(:,1);
mass1(:,i+1) = mass(:,i+1) - L1*q1plot(:,i+1);
mass2(:,i+1) = mass(:,i+1) - L2*q2plot(:,i+1);

R1 = sol(:,3:5);
R2 = sol(:,7:9);

detR1(i+1) = det(R1);
detR2(i+1) = det(R2);

rctime=controls(m1,m2,my,L1,L2,J1,J2,sol,d_trajectory(t));

q1dplot(:,i) = rctime(:,1);
q2dplot(:,i) = rctime(:,2);

errR1(i) = 0.5*trace(eye(3)-rctime(:,3:5)'*R1);
errR2(i) = 0.5*trace(eye(3)-rctime(:,6:8)'*R2);

% plot of the tracking problem
Q11(:,i+1) = mass1(:,i+1) + R1'*e1;
Q12(:,i+1) = mass1(:,i+1) - R1'*e1;
Q13(:,i+1) = mass1(:,i+1) + R1'*e2;
Q14(:,i+1) = mass1(:,i+1) - R1'*e2;

Q21(:,i+1) = mass2(:,i+1) + R2'*e1;
Q22(:,i+1) = mass2(:,i+1) - R2'*e1;
Q23(:,i+1) = mass2(:,i+1) + R2'*e2;
Q24(:,i+1) = mass2(:,i+1) - R2'*e2;

Xq1(:,i+1) = [Q11(1,i+1);Q13(1,i+1);Q12(1,i+1);Q14(1,i+1)];
Yq1(:,i+1) = [Q11(2,i+1);Q13(2,i+1);Q12(2,i+1);Q14(2,i+1)];
Zq1(:,i+1) = [Q11(3,i+1);Q13(3,i+1);Q12(3,i+1);Q14(3,i+1)];

Xq2(:,i+1) = [Q21(1,i+1);Q23(1,i+1);Q22(1,i+1);Q24(1,i+1)];
Yq2(:,i+1) = [Q21(2,i+1);Q23(2,i+1);Q22(2,i+1);Q24(2,i+1)];
Zq2(:,i+1) = [Q21(3,i+1);Q23(3,i+1);Q22(3,i+1);Q24(3,i+1)];

plot3(traj(1,:),traj(2,:),traj(3,:),'-')
axis([-3 3 -1 1 -9 0.5])
set(gca, 'zdir', 'reverse')
hold on;
plot3([mass(1,i+1),mass1(1,i+1)], [mass(2,i+1),mass1(2,i+1)], [mass(3,i+1),mass1(3,i+1)], 'k-');
axis([-3 3 -1 1 -9 0.5])
hold on;
plot3([mass(1,i+1),mass2(1,i+1)], [mass(2,i+1),mass2(2,i+1)], [mass(3,i+1),mass2(3,i+1)], 'k-');
axis([-3 3 -1 1 -9 0.5])
hold on;
plot3(mass1(1,i+1), mass1(2,i+1), mass1(3,i+1), 'd', 'MarkerSize',15, 'MarkerEdgeColor','k', 'MarkerFaceColor','k');
axis([-3 3 -1 1 -9 0.5])
hold on;
plot3(mass2(1,i+1), mass2(2,i+1), mass2(3,i+1), 'd', 'MarkerSize',15, 'MarkerEdgeColor','k', 'MarkerFaceColor','k');
axis([-3 3 -1 1 -9 0.5])
hold on;
plot3(mass(1,i+1), mass(2,i+1), mass(3,i+1), 'o', 'MarkerSize',8, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
axis([-3 3 -1 1 -9 0.5])

hold off;
pause(0.000000001)
