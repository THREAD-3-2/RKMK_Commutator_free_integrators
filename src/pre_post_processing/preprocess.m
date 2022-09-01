% Preprocessing : definition and inizialization of support variables to plot errors and trajectories

mass = zeros(3,N);
mass1 = mass;
mass2 = mass;

q1plot(:,1) = z0(:,11);
q2plot(:,1) = z0(:,13);
y1plot(:,1) = z0(1,1);
y2plot(:,1) = z0(2,1);
y3plot(:,1) = z0(3,1);
R1 = z0(:,3:5);
R2 = z0(:,7:9);

q1dplot = zeros(3,N-1);
q2dplot = zeros(3,N-1);

mass(:,1) = z0(:,1);
mass1(:,1) = z0(:,1) - L1*q10;
mass2(:,1) = z0(:,1) - L2*q20;
 
Q11 = zeros(3,N);
Q12 = zeros(3,N);
Q13 = zeros(3,N);
Q14 = zeros(3,N);
 
Q21 = Q11;
Q22 = Q11;
Q23 = Q11;
Q24 = Q11;
 
e1 = 0.15 * [1;0;0];
e2 = 0.15 * [0;1;0];
 
Q11(:,1) = mass1(:,1) + R1'*e1;
Q12(:,1) = mass1(:,1) - R1'*e1;
Q13(:,1) = mass1(:,1) + R1'*e2 ;
Q14(:,1) = mass1(:,1) - R1'*e2;
 
Q21(:,1) = mass2(:,1) + R2'*e1;
Q22(:,1) = mass2(:,1) - R2'*e1;
Q23(:,1) = mass2(:,1) + R2'*e2;
Q24(:,1) = mass2(:,1) - R2'*e2;
 
Xq1 = zeros(4,N);
Yq1 = Xq1;
Zq1 = Xq1;
Xq2 = Xq1;
Yq2 = Xq1;
Zq2 = Xq1;
 
Xq1(:,1) = [Q11(1,1);Q13(1,1);Q12(1,1);Q14(1,1)];
Yq1(:,1) = [Q11(2,1);Q13(2,1);Q12(2,1);Q14(2,1)];
Zq1(:,1) = [Q11(3,1);Q13(3,1);Q12(3,1);Q14(3,1)];
 
Xq2(:,1) = [Q21(1,1);Q23(1,1);Q22(1,1);Q24(1,1)];
Yq2(:,1) = [Q21(2,1);Q23(2,1);Q22(2,1);Q24(2,1)];
Zq2(:,1) = [Q21(3,1);Q23(3,1);Q22(3,1);Q24(3,1)];

traj = [time-2.5;zeros(size(time));((time-2.5).^2)-7];

figure;
plot3(traj(1,:),traj(2,:),traj(3,:),'-')
axis([-3 3 -1 1 -9 0.5])
set(gca, 'zdir', 'reverse')
hold on;

eynorm=zeros(1,N-1);

errR1=zeros(1,N-1);
errR2=zeros(1,N-1);

detR1=zeros(1,N);
detR2=zeros(1,N);
detR1(1)=det(R1);
detR2(1)=det(R2);
