%% »æÍ¼

close all
figure;
T = 6
t = (1:Time_sim/Tim_step)*Tim_step;
 plot(t,v0,'m--','linewidth',2);
hold on; plot(t, Velocity(:,1),'r','linewidth',2);
plot(t, Velocity(:,2),'b','linewidth',2);
plot(t, Velocity(:,3),'k','linewidth',2);
plot(t, Velocity(:,4),'g','linewidth',2);
plot(t, Velocity(:,5),'m','linewidth',2);
plot(t, Velocity(:,6),'r--','linewidth',2);
plot(t, Velocity(:,7),'b--','linewidth',2);
h = legend('0','1','2','3','4','5','6','7','8','Location','SouthEast');
set(h,'box','off'); box off;
xlabel('Time (s)');ylabel('Speed (m/s)');
axis([0 T 19 floor(max(max(Velocity))+1)])
set(gcf,'Position',[250 150 300 350]);

figure;
plot(t, Torque(:,1),'r','linewidth',2);hold on;
plot(t, Torque(:,2),'b','linewidth',2);hold on;
plot(t, Torque(:,3),'k','linewidth',2);hold on;
plot(t, Torque(:,4),'g','linewidth',2);hold on;
plot(t, Torque(:,5),'m','linewidth',2);hold on;
plot(t, Torque(:,6),'r--','linewidth',2);hold on;
plot(t, Torque(:,7),'b--','linewidth',2);hold on;
h = legend('1','2','3','4','5','6','7');
set(h,'box','off'); box off;
xlabel('Time (s)');ylabel('Torque (N)');
xlim([0 T])
set(gcf,'Position',[250 150 300 350]);

figure;
plot(t, (Eta*Torque(:,1)/R(1) - Ca(1)*Velocity(:,1).^2 - Mass(1)*g*f)/Mass(1),'r','linewidth',2);hold on;
plot(t, (Eta*Torque(:,2)/R(2) - Ca(2)*Velocity(:,2).^2 - Mass(2)*g*f)/Mass(2),'b','linewidth',2);hold on;
plot(t, (Eta*Torque(:,3)/R(3) - Ca(3)*Velocity(:,3).^2 - Mass(3)*g*f)/Mass(3),'k','linewidth',2);hold on;
plot(t, (Eta*Torque(:,4)/R(4) - Ca(4)*Velocity(:,4).^2 - Mass(4)*g*f)/Mass(4),'g','linewidth',2);hold on;
plot(t, (Eta*Torque(:,5)/R(5) - Ca(5)*Velocity(:,5).^2 - Mass(5)*g*f)/Mass(5),'m','linewidth',2);hold on;
plot(t, (Eta*Torque(:,6)/R(6) - Ca(6)*Velocity(:,6).^2 - Mass(6)*g*f)/Mass(6),'r--','linewidth',2);hold on;
plot(t, (Eta*Torque(:,7)/R(7) - Ca(7)*Velocity(:,7).^2 - Mass(7)*g*f)/Mass(7),'b--','linewidth',2);hold on;
h = legend('1','2','3','4','5','6','7','Location','NorthEast');
set(h,'box','off'); box off;
xlabel('Time (s)');ylabel('Acceleration (m/s^2)');
xlim([0 T])
set(gcf,'Position',[250 150 300 350]);

figure;
plot(t, x0 - Postion(:,1) - d,'r','linewidth',2);hold on;
plot(t, Postion(:,1) - d - Postion(:,2),'b','linewidth',2);hold on;
plot(t, Postion(:,2) - d - Postion(:,3),'k','linewidth',2);hold on;
plot(t, Postion(:,3) - d - Postion(:,4),'g','linewidth',2);hold on;
plot(t, Postion(:,4) - d - Postion(:,5),'m','linewidth',2);hold on;
plot(t, Postion(:,5) - d - Postion(:,6),'r--','linewidth',2);hold on;
plot(t, Postion(:,6) - d - Postion(:,7),'b--','linewidth',2);hold on;
h = legend('1','2','3','4','5','6','7');
set(h,'box','off'); box off;
xlabel('Time (s)');ylabel('Spacing error (m)');
xlim([0 T])
set(gcf,'Position',[250 150 300 350]);

figure;
plot(t, Postion(:,1)- (x0- d),'r','linewidth',2);hold on;
plot(t, Postion(:,2)-(x0- 2*d) ,'b','linewidth',2);hold on;
plot(t, Postion(:,3)-(x0- 3*d),'k','linewidth',2);hold on;
plot(t, Postion(:,4)-(x0- 4*d),'g','linewidth',2);hold on;
plot(t, Postion(:,5)-(x0- 5*d),'m','linewidth',2);hold on;
plot(t, Postion(:,6)-(x0- 6*d),'r--','linewidth',2);hold on;
plot(t, Postion(:,7)-(x0- 7*d),'b--','linewidth',2);hold on;
h = legend('1','2','3','4','5','6','7');
set(h,'box','off'); box off;
xlabel('Time (s)');ylabel('Spacing error (m)');
xlim([0 T])
set(gcf,'Position',[250 150 300 350]);