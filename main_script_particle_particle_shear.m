clear all
close all
cd 'C:\Users\chand\OneDrive - Wageningen University & Research\inter_particle_shear_github'
file_path='C:\Users\chand\OneDrive - UvA\paper_writing\inter_particle_shear\data_zenodo\data_set\Chandan-Particle_particle_shear_constant_w_mulitiple_level_06_01_2022\';
fig_save_loc='C:\Users\chand\OneDrive - UvA\paper_writing\inter_particle_shear\pictures\';
%Test 1
upper_particle_setup_dim=29.53/1000;
lower_particle_setup_dim=31.25/1000;
particle_dia_upper_para=18.16/1000;
particle_dia_upper_perp=18.55/1000;
particle_dia_lower_para=18.79/1000;
particle_dia_lower_perp=19.32/1000;
moment_arm=27/1000;

%upper particle
r_upper_mid=(particle_dia_upper_para+particle_dia_upper_perp)/2;
r_lower_mid=(particle_dia_lower_para+particle_dia_lower_perp)/2;
%650ml of milliQ water as medium

%% gap_1 %rep1
gap_no=1;
rep_no=1;
file_name='Test_59_9.txt';
file_loc=strcat(file_path,file_name);
cols_to_import=12;
Samplegap(gap_no).rep(rep_no).data=data_import(file_loc,cols_to_import);
Samplegap(gap_no).rep(rep_no).processed_data=process(Samplegap(gap_no).rep(rep_no).data,moment_arm,upper_particle_setup_dim,lower_particle_setup_dim,particle_dia_lower_perp,particle_dia_upper_perp);

close all
%% plot
f=figure
sf1=subplot(2,2,1)
for i=1:4
    name_tag=sprintf('%0.2e',Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean)
%     name_tag=strcat(num2str(Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean))
    plot(Samplegap(1).rep.processed_data.tofit(i, 1).x,Samplegap(1).rep.processed_data.tofit(i, 1).ver,'DisplayName',name_tag)
    hold on
end
yline(0,"DisplayName","Baseline")
hold off
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_v$ [N]",'Interpreter','latex')
lgd = legend('Location','south','Interpreter','latex');
lgd_title=sprintf('$v$ [m/s]')
lgd.Title.String = lgd_title;
lgd.FontSize = 10;
ax1 = gca;
ax1.LineWidth=1;
ax1.FontSize = 12;
ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')

sf2=subplot(2,2,2)
plot(Samplegap(1).rep.data(3).defangle.*moment_arm/1000,(Samplegap(1).rep.data(3).torque/moment_arm)/1000,'DisplayName','Clockwise')
% hold on
% plot(Samplegap(1).rep.data(5).defangle.*moment_arm/1000,(Samplegap(1).rep.data(5).torque/moment_arm)/1000,'DisplayName','Anti-Clockwise')
% hold off
[ymin,idx_min] = min((Samplegap(1).rep.data(3).torque/moment_arm)/1000) ;
[ymax,idx_max] = max((Samplegap(1).rep.data(3).torque/moment_arm)/1000) ; 
text(Samplegap(1).rep.data(3).defangle(idx_min).*moment_arm/1000,ymin,['t: ' num2str(round(ymin,4))]);
text(Samplegap(1).rep.data(3).defangle(idx_max).*moment_arm/1000,ymax,['p: ' num2str(round(ymax,4))]);
xlim([0.0045 0.0175])
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_h$ [N]",'Interpreter','latex')
ax2 = gca;
ax2.LineWidth=1;
ax2.FontSize = 12;
ax2.YAxis.Exponent = -3;
ax2.XAxis.Exponent = -3;
axis('square')


sf3=subplot(2,2,3)
plot(Samplegap(1).rep.data(3).defangle.*moment_arm/1000,(Samplegap(1).rep.data(3).torque/moment_arm)/1000,'DisplayName','C')
hold on
plot(Samplegap(1).rep.data(5).defangle.*moment_arm/1000,(Samplegap(1).rep.data(5).torque/moment_arm)/1000,'DisplayName','AC')
xline(Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,1)).*moment_arm/1000,'green','HandleVisibility','off')%'DisplayName',' Contact Clockwise')
xline(Samplegap(1).rep.data(5).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,2)).*moment_arm/1000,'m','HandleVisibility','off')%'DisplayName',' Contact Anti-Clockwise')
xline(mean([Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,1)),Samplegap(1).rep.data(5).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,2))]).*moment_arm/1000,'cyan','HandleVisibility','off')%'DisplayName','Mid Point')
yline(0,'HandleVisibility','off')%'DisplayName','Baseline')
hold off
xlim([0.0045 0.0175])
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_h$ [N]",'Interpreter','latex')
lgd = legend('Location','southwest');
lgd_title=sprintf('Dir')
lgd.Title.String = lgd_title;
lgd.FontSize = 10;
ax3 = gca;
ax3.LineWidth=1;
ax3.FontSize = 12;
ax3.YAxis.Exponent = -3;
ax3.XAxis.Exponent = -3;
axis('square')
axes('position',[.32 .33 .1 .1])
box on % put box around new pair of axes
plot(Samplegap(1).rep.data(3).defangle.*moment_arm/1000,(Samplegap(1).rep.data(3).torque/moment_arm)/1000,'DisplayName','C')
hold on
plot(Samplegap(1).rep.data(5).defangle.*moment_arm/1000,(Samplegap(1).rep.data(5).torque/moment_arm)/1000,'DisplayName','AC')
xline(mean([Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,1)),Samplegap(1).rep.data(5).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,2))]).*moment_arm/1000,'cyan','HandleVisibility','off')%'DisplayName','Mid Point')
yline(0,'HandleVisibility','off')%'DisplayName','Baseline')
hold off
ylim([-0.0003509122479578407 0.0003752653418182289])
xlim([0.010428698412224 0.010999628833344])
hold off
ax4 = gca;
ax4.LineWidth=1;
ax4.FontSize = 8;
ax4.YAxis.Exponent = -3;
ax4.XAxis.Exponent = -3;

sf4=subplot(2,2,4)
for i=1:4
    name_tag=sprintf('%0.2e',Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean)
%     name_tag=strcat(num2str(Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean))
    plot(Samplegap(1).rep.processed_data.tofit(i, 1).x,Samplegap(1).rep.processed_data.tofit(i, 1).F_hz,'DisplayName',name_tag)
    hold on
end
yline(0,"DisplayName","Baseline")
hold off
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_h$ [N]",'Interpreter','latex')
lgd = legend('Location','southwest','Interpreter','latex');
lgd_title=sprintf('$v$ [m/s]')
lgd.Title.String = lgd_title;
lgd.FontSize = 10;
ax5 = gca;
ax5.LineWidth=1;
ax5.FontSize = 12;
ax5.YAxis.Exponent = -3;
ax5.XAxis.Exponent = -3;
axis('square')

f.Position=[200 300 800 700]

AddLetters2Plots({ax1, ax2, ax3, ax5},{'a)', 'b)', 'c)', 'd)'}, 'HShift', [0.03, 0.03, 0.03, 0.03], 'VShift', 0)
% saveas(gcf,strcat(fig_save_loc,'exp_results_matlab.svg'))
% saveas(gcf,strcat(fig_save_loc,''exp_results_matlab.eps','epsc')

%% vertical forces for paper
%figure 8a)
f=figure
for i=1:4
    name_tag=sprintf('%0.2e',Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean)
%     name_tag=strcat(num2str(Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean))
    plot(Samplegap(1).rep.processed_data.tofit(i, 1).x,Samplegap(1).rep.processed_data.tofit(i, 1).ver,'DisplayName',name_tag)
    hold on
end
yline(0,"DisplayName","Baseline")
hold off
xlim([0 10e-3])
ylim([-10e-3 50e-3])
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_v$ [N]",'Interpreter','latex')
lgd = legend('Location','south','Interpreter','latex');
lgd_title=sprintf('v [m/s]')
lgd.Title.String = lgd_title;
lgd.FontSize = 18;
ax1 = gca;
ax1.LineWidth=2;
ax1.FontSize = 20;
ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')
f.Position=[200 300 560 500]
AddLetters2Plots({ax1},{'a)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,strcat(fig_save_loc,'ver_force_vs_disp.svg'))
saveas(gcf,strcat(fig_save_loc,'ver_force_vs_disp'),'epsc')
saveas(gcf,strcat(fig_save_loc,'ver_force_vs_disp.emf'))



%% horizontal forces for paper
% figure 1b) and 3a) raw
f=figure
plot(Samplegap(1).rep.data(3).defangle.*moment_arm/1000,(Samplegap(1).rep.data(3).torque/moment_arm)/1000,'DisplayName','Clockwise')
% hold on
% plot(Samplegap(1).rep.data(5).defangle.*moment_arm/1000,(Samplegap(1).rep.data(5).torque/moment_arm)/1000,'DisplayName','Anti-Clockwise')
% hold off
[ymin,idx_min] = min((Samplegap(1).rep.data(3).torque/moment_arm)/1000) ;
[ymax,idx_max] = max((Samplegap(1).rep.data(3).torque/moment_arm)/1000) ; 
text(Samplegap(1).rep.data(3).defangle(idx_min).*moment_arm/1000,ymin,['t: ' num2str(round(ymin,4))],'FontSize',18);
text(Samplegap(1).rep.data(3).defangle(idx_max).*moment_arm/1000,ymax,['p: ' num2str(round(ymax,4))],'FontSize',18);
xlim([0.0045 0.0175])
ylim([-0.0043 0.0043])
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_h$ [N]",'Interpreter','latex')
ax2 = gca;
ax2.LineWidth=2;
ax2.FontSize = 20;
ax2.YAxis.Exponent = -3;
ax2.XAxis.Exponent = -3;
axis('square')
f.Position=[200 300 560 500]
AddLetters2Plots({ax2},{'a)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,strcat(fig_save_loc,'hz_force_vs_disp_ampdiff.svg'))
saveas(gcf,strcat(fig_save_loc,'hz_force_vs_disp_ampdiff'),'epsc')
% eps file conversion after inkscape edit

%% contact lengths
% figure 4b)
f=figure
scatter(mean(abs(Samplegap(1).rep.processed_data.v_hz_mean),2),Samplegap(1).rep.processed_data.contact_def_len(:,1)*1000,80,'filled','DisplayName','$l_c$ C')
hold on
scatter(mean(abs(Samplegap(1).rep.processed_data.v_hz_mean),2),Samplegap(1).rep.processed_data.contact_def_len(:,2)*1000,80,'filled','DisplayName','$l_c$ AC')
scatter(mean(abs(Samplegap(1).rep.processed_data.v_hz_mean),2),Samplegap(1).rep.processed_data.lin_disp_m*1000,80,'d','filled','DisplayName','$l_g$')
hold off
xlabel("$v$ [m/s]",'Interpreter','latex')
ylabel("$l$ [mm]",'Interpreter','latex')
lgd = legend('Location','southeast');
set(lgd, 'Interpreter','latex')
lgd.FontSize = 18;
ax1 = gca;
ax1.LineWidth=2;
ax1.FontSize = 20;
% ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')
box on
f.Position=[200 300 560 500];
AddLetters2Plots({ax1},{'b)'}, 'HShift', 0.06, 'VShift', 0.1, 'FontSize', 18)
saveas(gcf,strcat(fig_save_loc,'contact_length_v_vel.svg'))
saveas(gcf,strcat(fig_save_loc,'contact_length_v_vel'),'epsc')
saveas(gcf,strcat(fig_save_loc,'contact_length_v_vel.emf'))

%% different velocities
% figure 4a)
f=figure
for i=1:4
    name_tag=sprintf('%0.2e',Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean)
%     name_tag=strcat(num2str(Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean))
    plot(Samplegap(1).rep.processed_data.tofit(i, 1).x,Samplegap(1).rep.processed_data.tofit(i, 1).F_hz,'DisplayName',name_tag,'LineWidth',2)
    hold on
end
yline(0,"DisplayName","Baseline")
rectangle('Position', [9.1e-3, -0.2e-3, 0.4e-3, 0.4e-3], 'EdgeColor', 'g');
hold off
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_h$ [N]",'Interpreter','latex')
lgd = legend('Location','southwest','Interpreter','latex');
lgd_title=sprintf('$v$ [m/s]')
lgd.Title.String = lgd_title;
lgd.FontSize = 18;
ax1 = gca;
ax1.LineWidth=2;
ax1.FontSize =20;
ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')
axes('position',[.6 .66 .2 .2])
box on % put box around new pair of axes
for i=1:4
    name_tag=sprintf('%0.2e',Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean)
%     name_tag=strcat(num2str(Samplegap(1).rep.processed_data.tofit(i, 1).v_hz_mean))
    plot(Samplegap(1).rep.processed_data.tofit(i, 1).x,Samplegap(1).rep.processed_data.tofit(i, 1).F_hz,'DisplayName',name_tag,'LineWidth',2)
    hold on
end
yline(0,'HandleVisibility','off')%'DisplayName','Baseline')
hold off
ylim([-0.2e-3 0.2e-3])
xlim([9.1e-3 9.5e-3])
hold off
ax4 = gca;
ax4.LineWidth=1;
ax4.FontSize = 14;
ax4.YAxis.Exponent = -3;
ax4.XAxis.Exponent = -3;
f.Position=[200 300 560 500];
AddLetters2Plots({ax1},{'a)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,strcat(fig_save_loc,'hz_force_diff_vel.svg'))
saveas(gcf,strcat(fig_save_loc,'hz_force_diff_vel'),'epsc')
saveas(gcf,strcat(fig_save_loc,'hz_force_diff_vel.emf'))

%% contact lengths both directions
% figure 3b)
f=figure
plot(Samplegap(1).rep.data(3).defangle.*moment_arm/1000,(Samplegap(1).rep.data(3).torque/moment_arm)/1000,'DisplayName','C')
hold on
plot(Samplegap(1).rep.data(5).defangle.*moment_arm/1000,(Samplegap(1).rep.data(5).torque/moment_arm)/1000,'DisplayName','AC')
% xline(Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,1)).*moment_arm/1000,'green','HandleVisibility','off')%'DisplayName',' Contact Clockwise')
% xline(Samplegap(1).rep.data(5).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,2)).*moment_arm/1000,'m','HandleVisibility','off')%'DisplayName',' Contact Anti-Clockwise')
xline(mean([Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,1)),Samplegap(1).rep.data(5).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,2))]).*moment_arm/1000,'black','HandleVisibility','off')%'DisplayName','Mid Point')
yline(0,'HandleVisibility','off')%'DisplayName','Baseline')
rectangle('Position', [0.010428698412224, -0.0003509122479578407, 0.010999628833344-0.010428698412224, 0.0003752653418182289-(-0.0003509122479578407)], 'EdgeColor', 'g');
hold off
xlim([0.0045 0.0175])
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_h$ [N]",'Interpreter','latex')
lgd = legend('Location','southwest');
lgd_title=sprintf('Dir')
lgd.Title.String = lgd_title;
lgd.FontSize = 18;
ax3 = gca;
ax3.LineWidth=2;
ax3.FontSize = 20;
ax3.YAxis.Exponent = -3;
ax3.XAxis.Exponent = -3;
axis('square')
axes('position',[.6 .66 .2 .2])
box on % put box around new pair of axes
plot(Samplegap(1).rep.data(3).defangle.*moment_arm/1000,(Samplegap(1).rep.data(3).torque/moment_arm)/1000,'DisplayName','C')
hold on
plot(Samplegap(1).rep.data(5).defangle.*moment_arm/1000,(Samplegap(1).rep.data(5).torque/moment_arm)/1000,'DisplayName','AC')
xline(mean([Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,1)),Samplegap(1).rep.data(5).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,2))]).*moment_arm/1000,'black','HandleVisibility','off')%'DisplayName','Mid Point')
yline(0,'HandleVisibility','off')%'DisplayName','Baseline')
hold off
ylim([-0.0003509122479578407 0.0003752653418182289])
xlim([0.010428698412224 0.010999628833344])
hold off
ax4 = gca;
ax4.LineWidth=1;
ax4.FontSize = 14;
ax4.YAxis.Exponent = -3;
ax4.XAxis.Exponent = -3;
f.Position=[200 300 560 500]
AddLetters2Plots({ax3},{'b)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
% saveas(gcf,'M:\paper_writing\inter_particle_shear\hz_force_vs_disp_diff_direction.svg')
% saveas(gcf,'M:\paper_writing\inter_particle_shear\hz_force_vs_disp_diff_direction','epsc')


%% 
figure
plot(Samplegap(1).rep.data(3).defangle.*moment_arm/1000,(Samplegap(1).rep.data(3).torque/moment_arm)/1000,'DisplayName','Clockwise')
hold on
% plot(Samplegap.rep.data(5).defangle.*moment_arm/1000,(Samplegap.rep.data(5).torque/moment_arm)/1000,'DisplayName','Anti-Clockwise')
xline(Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,1)).*moment_arm/1000,'green','DisplayName',' Contact Clockwise')
xline(Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_stop_index(1,1)).*moment_arm/1000,'m','DisplayName',' Contact lost')
xline(mean([Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_start_index(1,1)),Samplegap(1).rep.data(3).defangle(Samplegap(1).rep.processed_data.filter.lim_stop_index(1,1))]).*moment_arm/1000,'cyan','DisplayName','Mid Point')
yline(0,'DisplayName','Baseline')
hold off
xlim([0.0045 0.0175])
xlabel('Displacement(m)')
ylabel('Raw Horizontal Force(N)')
lgd = legend('Location','south');
lgd.Title.String = 'Rotational direction';


figure
plot(Samplegap(1).rep.processed_data.tofit(1, 1).x/(mean([particle_dia_lower_perp particle_dia_upper_perp])/2),Samplegap(1).rep.processed_data.tofit(1, 1).F_hz/Samplegap(1).rep.processed_data.overlap(1,1))
hold on
xline((Samplegap(1).rep.processed_data.lin_disp_m(1)/2)/(mean([particle_dia_lower_perp particle_dia_upper_perp])/2))
yline(0)
hold off
xlabel('^{Displacement}/_{Avg. Radius}')
ylabel('^{Horizontal Force(N)}/_{Max. overlap(m)}')


% horizontal force vs displacement
f=figure
plot(Samplegap(1).rep.processed_data.tofit(1, 1).x,Samplegap(1).rep.processed_data.tofit(1, 1).F_hz,'LineWidth',2)
hold on
% xline((Samplegap.rep.processed_data.lin_disp_m(1)/2)/(mean([particle_dia_lower_perp particle_dia_upper_perp])/2))
yline(0)
hold off
box on
xlabel('$\theta$ [m]','Interpreter','latex')
ylabel('$F_h$ [N]','Interpreter','latex')
ax = gca;
ax.LineWidth=2;
ax.FontSize = 25;
ax.YAxis.Exponent = -3;
ax.XAxis.Exponent = -3;
axis('square')
f.Position=[200 300 560 500]
% hFig = figure(1);
% set(gcf,'PaperPosition',[10 10 20 20])
saveas(gcf,strcat(fig_save_loc,'Hz_force_vs_disp.svg'))


% vertical force vs displacement
f=figure
plot(Samplegap(1).rep.processed_data.tofit(1, 1).x,Samplegap(1).rep.processed_data.tofit(1, 1).ver,'LineWidth',2)
hold on
% xline(Samplegap.rep.processed_data.lin_disp_m(1)/2)
yline(0)
hold off
box on
xlabel('$\theta$ [m]','Interpreter','latex')
ylabel('$F_v$ [N]','Interpreter','latex')
ax = gca;
ax.LineWidth=2;
ax.FontSize = 25;
ax.YAxis.Exponent = -3;
ax.XAxis.Exponent = -3;
axis('square')
f.Position=[200 300 560 500]
saveas(gcf,strcat(fig_save_loc,'Ver_force_vs_disp.svg'))

% comparing hz force signals

figure
vel_mean=mean(abs(Samplegap(1).rep.processed_data.v_hz_mean),2);
seperation_mean=mean(abs(Samplegap(1).rep.processed_data.sep_clock),2);
scatter(vel_mean,seperation_mean,'filled')
xlabel('Rotational velocity(m/s)')
ylabel('Difference detachment and contact(m)')

% change in contact length with velocity
figure
vel_mean=mean(abs(Samplegap(1).rep.processed_data.v_hz_mean),2);
seperation_mean_ratio=mean(abs(Samplegap(1).rep.processed_data.sep_clock),2)./Samplegap(1).rep.processed_data.lin_disp_m;
scatter(vel_mean,seperation_mean,'filled')
xlabel('Rotational velocity(m/s)')
ylabel('Change in contact length(m)/geometric contact length')


figure
scatter(vel_mean,Samplegap(1).rep.processed_data.avg_contact_len'./Samplegap(1).rep.processed_data.lin_disp_m)
xlabel('Rotational velocity(m/s)')
ylabel('$\frac{force contact length}{geometric contact length}$','Interpreter','latex','FontSize',20)

% data export to fit in x/x_e vs v/v_char
lin_disp_m=Samplegap.rep(1).processed_data.lin_disp_m;
save('vel_variables_1.mat','seperation_mean','vel_mean', 'lin_disp_m')

lin_disp_m_array(gap_no,rep_no)=mean(lin_disp_m);
lin_disp_m_array_error(gap_no,rep_no)=max(lin_disp_m)-min(lin_disp_m);



figure
peak_first=mean(Samplegap.rep.processed_data.max_amp,2);
peak_second=abs(mean(Samplegap.rep.processed_data.min_amp,2));
scatter(Samplegap.rep.processed_data.v_hz_mean(:,1),peak_first,'filled','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'DisplayName','First Peak')
hold on
scatter(Samplegap.rep.processed_data.v_hz_mean(:,2),peak_second,'filled','MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'DisplayName','Second Peak')
hold off
xlabel('Rotational velocity(m/s)')
ylabel('Amplitudes(N)')


peak_first=mean(Samplegap(1).rep.processed_data.max_amp,2);
peak_second=(mean(Samplegap(1).rep.processed_data.min_amp,2));
figure
scatter(Samplegap(1).rep.processed_data.v_hz_mean(:,1),peak_first+peak_second,'filled','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'DisplayName','First Peak')

xlabel('Rotational velocity(m/s)')
ylabel('Diff. in Amplitudes(N)')


figure
mean_geom_cycle_center=mean(Samplegap(1).rep.processed_data.lin_disp_m)/2;
x_zero_mean=mean(Samplegap(1).rep.processed_data.x_zero_exp,2)
scatter(Samplegap(1).rep.processed_data.v_hz_mean(:,1),abs(mean_geom_cycle_center-x_zero_mean),'filled')%,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410])

xlabel('Rotational velocity(m/s)')
ylabel('Central Shift(m)')

%% horizontal force signal analysis
figure
peak_first=mean(Samplegap(1).rep.processed_data.max_amp,2);
peak_second=(mean(Samplegap(1).rep.processed_data.min_amp,2));
subplot(1,2,1)
scatter(Samplegap(1).rep.processed_data.v_hz_mean(:,1),peak_first+peak_second,'filled','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410])
xlabel('Rotational velocity(m/s)')
ylabel('Diff. in Amplitudes(N)')

subplot(1,2,2)
scatter(vel_mean,Samplegap(1).rep.processed_data.avg_contact_len'./Samplegap(1).rep.processed_data.lin_disp_m,'filled','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410])
xlabel('Rotational velocity(m/s)')
ylabel('$\frac{force contact length}{geometric contact length}$','Interpreter','latex','FontSize',20)


%% horizontal force signal analysis
figure
peak_first=mean(Samplegap(1).rep.processed_data.max_amp,2);
peak_second=(mean(Samplegap(1).rep.processed_data.min_amp,2));
subplot(1,3,1)
scatter(Samplegap(1).rep.processed_data.v_hz_mean(:,1),peak_first+peak_second,'filled','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410])

xlabel('Rotational velocity(m/s)')
ylabel('Difference in Amplitudes(N)')

subplot(1,3,2)
scatter(vel_mean,Samplegap(1).rep.processed_data.avg_contact_len','filled','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410])
xlabel('Rotational velocity(m/s)')
ylabel('Force contact length (m)')

subplot(1,3,3)
for i=1:1
scatter(Samplegap(i).rep.processed_data.v_hz_mean(:,1),Samplegap(i).rep.processed_data.x_zero_exp(:,1),'filled','MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'DisplayName','Clockwise');
hold on;
scatter(Samplegap(i).rep.processed_data.v_hz_mean(:,end),Samplegap(i).rep.processed_data.x_zero_exp(:,end),'MarkerEdgeColor',[0 0.4470 0.7410],'DisplayName','Anticlockwise');
end
hold off
xlabel('Rotational velocity(m/s)')
ylabel('Asymmetry in zero crossing (m)')
legend

%% for ovelap fit
disp_mean(1)=mean(Samplegap(gap_no).rep.processed_data.lin_disp_m);
gap_mean(1)=mean(Samplegap(gap_no).rep.data(1).gap);

%% gap_2 %rep1
gap_no=2;
rep_no=1;
file_name='Test_59_8.txt';
file_loc=strcat(file_path,file_name);
cols_to_import=12;
Samplegap(gap_no).rep(rep_no).data=data_import(file_loc,cols_to_import);

Samplegap(gap_no).rep(rep_no).processed_data=process(Samplegap(gap_no).rep(rep_no).data,moment_arm,upper_particle_setup_dim,lower_particle_setup_dim,particle_dia_lower_perp,particle_dia_upper_perp);


%% for ovelap fit
disp_mean(2)=mean(Samplegap(gap_no).rep.processed_data.lin_disp_m);
gap_mean(2)=mean(Samplegap(gap_no).rep.data(1).gap);


para_overlap.gap1_2=over_lap(disp_mean(1),disp_mean(2),particle_dia_upper_para,particle_dia_lower_para,gap_mean(1),gap_mean(2))

%% gap_3 %rep1
gap_no=3;
rep_no=1;
file_name='Test_59_7.txt';
file_loc=strcat(file_path,file_name);
cols_to_import=12;
Samplegap(gap_no).rep(rep_no).data=data_import(file_loc,cols_to_import);
Samplegap(gap_no).rep(rep_no).processed_data=process(Samplegap(gap_no).rep(rep_no).data,moment_arm,upper_particle_setup_dim,lower_particle_setup_dim,particle_dia_lower_perp,particle_dia_upper_perp);

%% for ovelap fit
disp_mean(3)=mean(Samplegap(gap_no).rep.processed_data.lin_disp_m);
gap_mean(3)=mean(Samplegap(gap_no).rep.data(1).gap);


para_overlap.gap2_3=over_lap(disp_mean(2),disp_mean(3),particle_dia_upper_para,particle_dia_lower_para,gap_mean(2),gap_mean(3))


%% gap_4 %rep1
gap_no=4;
rep_no=1;
file_name='Test_59_6.txt';
file_loc=strcat(file_path,file_name);
cols_to_import=12;
Samplegap(gap_no).rep(rep_no).data=data_import(file_loc,cols_to_import);
Samplegap(gap_no).rep(rep_no).processed_data=process(Samplegap(gap_no).rep(rep_no).data,moment_arm,upper_particle_setup_dim,lower_particle_setup_dim,particle_dia_lower_perp,particle_dia_upper_perp);


% data export to fit in x/x_e vs v/v_char
lin_disp_m=Samplegap(gap_no).rep(1).processed_data.lin_disp_m;
save('vel_variables_4.mat','seperation_mean','vel_mean', 'lin_disp_m')

lin_disp_m_array(gap_no,rep_no)=mean(lin_disp_m);
lin_disp_m_array_error(gap_no,rep_no)=max(lin_disp_m)-min(lin_disp_m);

%% gap_5 %rep1
gap_no=5;
rep_no=1;
file_name='Test_59_5.txt';
file_loc=strcat(file_path,file_name);
cols_to_import=12;
Samplegap(gap_no).rep(rep_no).data=data_import(file_loc,cols_to_import);
Samplegap(gap_no).rep(rep_no).processed_data=process(Samplegap(gap_no).rep(rep_no).data,moment_arm,upper_particle_setup_dim,lower_particle_setup_dim,particle_dia_lower_perp,particle_dia_upper_perp);

% data export to fit in x/x_e vs v/v_char
lin_disp_m=Samplegap(gap_no).rep(1).processed_data.lin_disp_m;
save('vel_variables_5.mat','seperation_mean','vel_mean', 'lin_disp_m')

lin_disp_m_array(gap_no,rep_no)=mean(lin_disp_m);
lin_disp_m_array_error(gap_no,rep_no)=max(lin_disp_m)-min(lin_disp_m);



%% for paper
% figure 7a)
f=figure
gap_no=5;
for i=1:gap_no
    over_lap_mean=(mean(mean(Samplegap(i).rep.processed_data.overlap)))-1e-4;
    name_tag=sprintf('%0.2e',over_lap_mean);
    start_angle=(Samplegap(i).rep.processed_data.def_angle_start(1,1)*moment_arm/1000);
    plot(Samplegap(i).rep.processed_data.tofit(1, 1).x+start_angle,Samplegap(i).rep.processed_data.tofit(1, 1).F_hz,'DisplayName',name_tag,'LineWidth',2)
    hold on
end
yline(0,'HandleVisibility','off')
hold off
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_h$ [N]",'Interpreter','latex')
lgd = legend('Location','southwest','Interpreter','latex');
lgd_title='$\delta_{max}$ [m]';
lgd.Title.String = lgd_title;
lgd.FontSize = 17;
ax1 = gca;
ax1.LineWidth=2;
ax1.FontSize =20;
ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')
f.Position=[200 300 560 500];
AddLetters2Plots({ax1},{'a)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,strcat(fig_save_loc,'hz_force_diff_overlap.svg'))
saveas(gcf,strcat(fig_save_loc,'hz_force_diff_overlap'),'epsc')
saveas(gcf,strcat(fig_save_loc,'hz_force_diff_overlap.emf'))

%figure 7b)
f=figure
gap_no=5;
for i=1:gap_no
    over_lap_mean=(mean(mean(Samplegap(i).rep.processed_data.overlap)))-1e-4;
     name_tag=sprintf('%0.2e',over_lap_mean);
    scatter(Samplegap(i).rep(1).processed_data.v_hz_mean_avg,Samplegap(i).rep(1).processed_data.contact_def_len_mean/Samplegap(i).rep(1).processed_data.contact_def_len_geometric,80,'filled','DisplayName',name_tag)
    hold on
end
hold off
ylim([0.94 1])
ytickformat('%.2f')
xlabel("$v$ [m/s]",'Interpreter','latex')
ylabel("${l_c}/{l_g}$",'Interpreter','latex')
lgd = legend('Location','southeast','Interpreter','latex');
lgd_title='$\delta_{max} [m]$';
lgd.Title.String = lgd_title;
lgd.FontSize = 18;
ax1 = gca;
ax1.LineWidth=2;
ax1.FontSize = 20;
% ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')
box on
f.Position=[200 300 560 500]
AddLetters2Plots({ax1},{'b)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,strcat(fig_save_loc,'l_c_by_lg.svg'))
saveas(gcf,strcat(fig_save_loc,'l_c_by_lg'),'epsc')
saveas(gcf,strcat(fig_save_loc,'l_c_by_lg.emf'))

% figure 8c)
f=figure
for i=1:gap_no
    over_lap_mean=(mean(mean(Samplegap(i).rep.processed_data.overlap)))-1e-4;
    name_tag=sprintf('%0.2e',over_lap_mean);
    scatter(Samplegap(i).rep(1).processed_data.v_hz_mean_avg,Samplegap(i).rep.processed_data.mu_exp_mean,80,'filled','DisplayName',name_tag)
    hold on
end
hold off
box on
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel("$v$ [m/s]",'Interpreter','latex')
ylabel("$\mu$",'Interpreter','latex')
lgd = legend('Location','northeast','Interpreter','latex');
lgd_title='$\delta_{max} [m]$';
lgd.Title.String = lgd_title;
lgd.FontSize = 18;
ax1 = gca;
ax1.LineWidth=2;
ax1.FontSize = 20;
ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')
box on
f.Position=[200 300 560 500]
AddLetters2Plots({ax1},{'c)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,strcat(fig_save_loc,'mu_vs_v.svg'))
saveas(gcf,strcat(fig_save_loc,'mu_vs_v'),'epsc')
saveas(gcf,strcat(fig_save_loc,'mu_vs_v.emf'))

