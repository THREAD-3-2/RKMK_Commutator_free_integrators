% plots

figure;
subplot(2,1,1)
plot(timeSol, 1-detR1ref, 'b-', timeSol, 1-detR2ref, 'k*', 'linewidth', 1.5);
xlabel('Time', 'fontsize', 20)
ylabel('Determinant', 'fontsize', 20)
legend('R_1', 'R_2')
title('Ode45', 'fontsize', 15)
subplot(2,1,2)
plot(time, 1-detR1, 'b-', time, 1-detR2, 'k*', 'linewidth', 1.5 );
xlabel('Time', 'fontsize', 20)
ylabel('Determinant', 'fontsize', 20)
legend('R_1', 'R_2')
title('RKMK 4', 'fontsize', 15)
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
h=sgtitle('$1-detR_i(t)$', 'interpreter', 'Latex');
h.FontSize = 20;


figure;
subplot(2,1,1)
plot(timeSol, 1-vecnorm(q1ref).^2, 'b-', timeSol, 1-vecnorm(q2ref).^2, 'k-', 'linewidth', 1.5);
xlabel('Time', 'fontsize', 20)
ylabel('Norm', 'fontsize', 20)
h=legend('q_1', 'q_2');
h.FontSize=15;
title('Ode45', 'fontsize', 15)
subplot(2,1,2)
plot(time, 1-vecnorm(q1plot).^2, 'b-',time, 1-vecnorm(q2plot).^2, 'k-', 'linewidth', 1.5 );
xlabel('Time', 'fontsize', 20)
ylabel('Norm', 'fontsize', 20)
h=legend('q_1', 'q_2');
h.FontSize=15;
title('RKMK 4', 'fontsize', 15)
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
h=sgtitle('$1-q_i^T(t)q_i(t)$', 'interpreter', 'Latex');
h.FontSize = 20;


figure;
plot(time(2:end), eynorm, 'k-', 'linewidth', 1.5);
xlabel('Time', 'fontsize', 20)
h=ylabel('$\|y(t)-y_d(t)\|_2$', 'interpreter', 'Latex');
h.FontSize = 20;
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;


figure;
subplot(3,1,1)
plot(time, y1plot, 'b-', time, traj(1,:), 'r-', 'linewidth', 1.5);
h=legend('$y^{(1)}$', '$y_d^{(1)}$', 'interpreter', 'Latex', 'location', 'southeast');
h.FontSize=15;
title('Component 1', 'fontsize', 15);
% h.FontSize=20;
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
subplot(3,1,2)
plot(time, y2plot, 'b-', time, traj(2,:), 'r-', 'linewidth', 1.5);
h=legend('$y^{(2)}$', '$y_d^{(2)}$', 'interpreter', 'Latex', 'location', 'northeast');
h.FontSize=15;
title('Component 2', 'fontsize', 15);
% h.FontSize=15;
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
subplot(3,1,3)
plot(time, y3plot, 'b-', time, traj(3,:), 'r-', 'linewidth', 1.5);
set(gca, 'ydir', 'reverse')
h=legend('$y^{(3)}$', '$y_d^{(3)}$', 'interpreter', 'Latex', 'location', 'northeast');
h.FontSize=15;
xlabel('Time', 'fontsize', 20);
title('Component 3', 'fontsize', 15);
% h.FontSize=15;
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;


figure;
range = [5 500 1200 2000];
plot3(traj(1,:),traj(2,:),traj(3,:),'-')
axis([-3 3 -1 1 -9 0.5])
set(gca, 'zdir', 'reverse')
axis off
hold on;
quiver3([0;0;0],[0;0;0],[-1;-1;-1],0.8*[1;0;0],0.5*[0;1;0],[0;0;1]);
hold on
for i = range
    plot3([mass(1,i),mass1(1,i)], [mass(2,i),mass1(2,i)], [mass(3,i),mass1(3,i)], 'k-');
    axis([-3 3 -1 1 -9 0.5])
    axis off
    hold on;
    plot3([mass(1,i),mass2(1,i)], [mass(2,i),mass2(2,i)], [mass(3,i),mass2(3,i)], 'k-');
    axis([-3 3 -1 1 -9 0.5])
    axis off
    hold on;
    plot3(mass1(1,i), mass1(2,i), mass1(3,i), 'd', 'MarkerSize',15, 'MarkerEdgeColor','k', 'MarkerFaceColor','k');
    axis([-3 3 -1 1 -9 0.5])
    axis off
    hold on;
    plot3(mass2(1,i), mass2(2,i), mass2(3,i), 'd', 'MarkerSize',15, 'MarkerEdgeColor','k', 'MarkerFaceColor','k');
    axis([-3 3 -1 1 -9 0.5])
    axis off
    hold on;
    plot3(mass(1,i), mass(2,i), mass(3,i), 'o', 'MarkerSize',8, 'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    axis([-3 3 -1 1 -9 0.5])
    axis off
    hold on;
end
