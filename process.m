function processed_data=process(data,moment_arm,upper_particle_setup_dim,lower_particle_setup_dim,particle_dia_lower_perp,particle_dia_upper_perp)

start_interval=3;
step_by=4;
reversal_int=2;
s=size(data);
cou_nt=1;
num_speeds=4;
col=hsv(num_speeds);
window=3;
figure(1)
for i=start_interval:step_by:s(2)
    %clockwise % 3rd minesion 1 is clockwise
    F_tan(cou_nt,:,1)=(data(i).torque/moment_arm)/1000; %tangentaial force in N
    F_norm(cou_nt,:,1)=data(i).normal_force;
    direction(cou_nt,1)=0;
    time_exp(cou_nt,:,1)=data(i).interval_time;
    def_ang(cou_nt,:,1)=data(i).defangle;
    v_hz(cou_nt,:,1)=(data(i).rot_speed/1000)*moment_arm; % linear hz velocity in m/s
    overlap_avg(cou_nt,1)=upper_particle_setup_dim+lower_particle_setup_dim-(mean(data(i).gap)/1000);
    dtorque_dtheta(cou_nt,:,1)=diff(squeeze(F_tan(cou_nt,:,1)))./diff(squeeze(def_ang(cou_nt,:,1))); %N/mrad
    dtorque_dtheta_movmean(cou_nt,:,1)=movmean(dtorque_dtheta(cou_nt,:,1),window);
    plot (squeeze(def_ang(cou_nt,:,1)),squeeze(F_tan(cou_nt,:,1)),"DisplayName",num2str(cou_nt),"Color",col(cou_nt,:),"LineStyle","-")
    hold on
    plot (squeeze(def_ang(cou_nt,:,1)),squeeze(F_norm(cou_nt,:,1)),"DisplayName",num2str(cou_nt),"Color",col(cou_nt,:),"LineStyle","-")
    
    %anticlock % 3rd minesion 2 is anticlockwise
    F_tan(cou_nt,:,2)=(data(i+reversal_int).torque/moment_arm)/1000; %tangentaial force in N
    F_norm(cou_nt,:,2)=data(i+reversal_int).normal_force;
    direction(cou_nt,:,2)=1;
    time_exp(cou_nt,:,2)=data(i+reversal_int).interval_time;
    def_ang(cou_nt,:,2)=data(i+reversal_int).defangle;
    v_hz(cou_nt,:,2)=(data(i+reversal_int).rot_speed/1000)*moment_arm; %linear hz velocity in m/s
    overlap_avg(cou_nt,2)=upper_particle_setup_dim+lower_particle_setup_dim-(mean(data(i+reversal_int).gap)/1000);
    dtorque_dtheta(cou_nt,:,2)=diff(squeeze(F_tan(cou_nt,:,2)))./diff(squeeze(def_ang(cou_nt,:,2)));
    dtorque_dtheta_movmean(cou_nt,:,2)=movmean(dtorque_dtheta(cou_nt,:,2),window);
    plot (squeeze(def_ang(cou_nt,:,2)),squeeze(F_tan(cou_nt,:,2)),"DisplayName",num2str(cou_nt),"Color",col(cou_nt,:),"LineStyle","-.")
    plot (squeeze(def_ang(cou_nt,:,2)),squeeze(F_norm(cou_nt,:,2)),"DisplayName",num2str(cou_nt),"Color",col(cou_nt,:),"LineStyle","-.")
    
    cou_nt=cou_nt+1;
end
yline(0)
hold off
xlabel("Angular Displacement(mrad)")
ylabel("Forces in Hz direction")

figure
for i=1:cou_nt-1
    plot ((def_ang(i,:,1))*moment_arm/1000,F_tan(i,:,1),"DisplayName",num2str(cou_nt),"Color",col(i,:),"LineStyle","-")
    hold on
    plot ((def_ang(i,:,2))*moment_arm/1000,F_tan(i,:,2),"DisplayName",num2str(cou_nt),"Color",col(i,:),"LineStyle","-.")
end
yline(0)
hold off
xlabel("Displacement(m)")
ylabel("Forces in Hz direction(N)")

%% histogram of Tangential force data
window_hist=100;

%clockwise
% for j=1:2
% figure
% for i=1:(cou_nt-1)
%     subplot(2,2,i)
%     hist(squeeze(F_tan(i,:,1)),window_hist)
%     xlabel("Tangential Force (N)")
%     ylabel("frequency")
% end
%     if j==1 
%         direc="Clockwise";
%     else
%         direc="Anticlockwise"; 
%     end
% sgtitle(strcat("Tangential Force Distribution",direc))
% end

%binning

% for j=1:2
%     for i=1:(cou_nt-1)
%         [Fbin(i).N(j,:),Fbin(i).edges(j,:)] = histcounts(F_tan(i,:,1),window_hist);
%         [mode,I]=max(Fbin(i).N(j,:));
%         int_start_Tor(i,j)=Fbin(i).edges(1,I);
%         int_prev_Tor(i,j)=Fbin(i).edges(1,I-1);
%         int_end_Tor(i,j)=Fbin(i).edges(1,I+1);
%     end
% end

%% correcting for slope in baseline

% i=1; %speed level 1
% j=2; %clockwise

    for i=1:(cou_nt-1)
        % torque
%         x_data=[def_ang(i,100:400,j) def_ang(i,end-400:end,j)-def_ang(i,end-400,j)+def_ang(i,400,j)];
%         y_data=[F_tan(i,100:400,j) F_tan(i,end-400:end,j)];
%         fit_baseline_tor=polyfit(x_data,y_data,1);
%         figure
%         scatter(x_data,y_data)
%         hold on
%         plot(x_data,polyval(fit_baseline_tor,x_data))
%         hold off
%         title('Torques')
%         clear  x_data y_data
%         figure
%         plot(squeeze(F_tan(i,:,j)))
%         hold on
%         F_tan_c(i,:,j)=F_tan(i,:,j)-polyval(fit_baseline_tor,def_ang(i,:,j));
%         plot(squeeze(F_tan_c(i,:,j)))
%         hold off
        % vertical forces
%         x_data=[def_ang(i,100:400,j) def_ang(i,end-400:end,j)-def_ang(i,end-400,j)+def_ang(i,400,j)];
%         y_data=[F_norm(i,100:400,j) F_norm(i,end-400:end,j)];
%         x_data=[def_ang(i,end-400:end,j)-def_ang(i,end-400,j)];
%         y_data=[F_norm(i,end-400:end,j)];
        j=1;
        x_data=[def_ang(i,100:400,j) def_ang(i,end-700:end-100,j)];
        y_data=[F_norm(i,100:400,j) F_norm(i,end-700:end-100,j)];
        fit_baseline_norm=polyfit(x_data,y_data,1);
        
%         figure
%         scatter(x_data,y_data)
%         hold on
%         plot(x_data,polyval(fit_baseline_norm,x_data))
%         hold off
%         title('Vertical Forces')
        
        clear x_data y_data
        
%         figure
%         plot(squeeze(F_norm(i,:,j)))
%         hold on
%         yline(0)
         F_norm_c(i,:,j)=F_norm(i,:,j)-polyval(fit_baseline_norm,def_ang(i,:,j));
%         plot(squeeze(F_norm_c(i,:,j)))
%         hold off
        
        F_norm_o(i,:,j)=F_norm(i,:,j);
        F_norm(i,:,j)=F_norm_c(i,:,j);
        F_norm_max(i,j)=max(F_norm(i,:,j));
        j=2;
        
        x_data=[def_ang(i,100:700,j) def_ang(i,end-400:end-100,j)];
        y_data=[F_norm(i,100:700,j) F_norm(i,end-400:end-100,j)];
        fit_baseline_norm=polyfit(x_data,y_data,1);
        
%         figure
%         scatter(x_data,y_data)
%         hold on
%         plot(x_data,polyval(fit_baseline_norm,x_data))
%         hold off
%         title('Vertical Forces')

        clear x_data y_data
        
%         figure
%         plot(squeeze(F_norm(i,:,j)))
%         hold on
%         yline(0)
         F_norm_c(i,:,j)=F_norm(i,:,j)-polyval(fit_baseline_norm,def_ang(i,:,j));
%         plot(squeeze(F_norm_c(i,:,j)))
%         hold off

        F_norm_o(i,:,j)=F_norm(i,:,j);
        F_norm(i,:,j)=F_norm_c(i,:,j);
        F_norm_max(i,j)=max(F_norm(i,:,j));
end

figure
j=1
for i=1:cou_nt-1
    scatter(def_ang(i,:,j),F_norm(i,:,j),'filled','MarkerFaceColor',col(i,:))
    hold on
end
yline(0)
j=2
for i=1:cou_nt-1
    scatter(def_ang(i,:,j),F_norm(i,:,j),'d','filled','MarkerFaceColor',col(i,:))
    hold on
end
% yline(0)
hold off
xlabel("$\theta$ [mrad]",'Interpreter','latex')
ylabel("$F_v$ [N]",'Interpreter','latex')

%% std deviation filtering

j=1
    for i=1:(cou_nt-1)
%         figure
        baseline_avg(i,j)=mean(squeeze(F_tan(i,end-700:end-100,j)));
        std_dev(i,j)=std(squeeze(F_tan(i,end-700:end-100,j)));
%         F_tan(1,:,j)=F_tan(1,:,j)-baseline_avg(i,j);
%         baseline_avg(i,j)=mean(squeeze(F_tan(i,end-700:end-100,j)));
        baseline_avg_vert(i,j)=mean(squeeze(F_norm(i,end-700:end-100,j)));
        std_dev_vert(i,j)=std(squeeze(F_norm(i,end-700:end-100,j)));
        
%         plot(squeeze(F_tan(i,100:400,j)))
%         yline(baseline_avg(i,j))
%         yline(baseline_avg(i,j)+2*std_dev(i,j),'-.')
%         yline(baseline_avg(i,j)-2*std_dev(i,j),'--')

        
    end
    
j=2
    for i=1:(cou_nt-1)
%         figure
        baseline_avg(i,j)=mean(squeeze(F_tan(i,100:700,j)));
        std_dev(i,j)=std(squeeze(F_tan(i,100:700,j)));
%         F_tan(1,:,j)=F_tan(1,:,j)-baseline_avg(i,j);
%         baseline_avg(i,j)=mean(squeeze(F_tan(i,100:700,j)));
        baseline_avg_vert(i,j)=mean(squeeze(F_norm(i,100:700,j)));
        std_dev_vert(i,j)=std(squeeze(F_norm(i,100:700,j)));
        
%         plot(squeeze(F_tan(i,100:400,j)))
%         yline(baseline_avg(i,j))
%         yline(baseline_avg(i,j)+2*std_dev(i,j),'-.')
%         yline(baseline_avg(i,j)-2*std_dev(i,j),'--')
    end
    
% mean_baseline_F_tan=mean(baseline_avg);
% F_tan(:,:,1)=F_tan(:,:,1)-mean_baseline_F_tan(1);
% F_tan(:,:,2)=F_tan(:,:,2)-mean_baseline_F_tan(2);
% 
% 
% j=1;
% for i=1:(cou_nt-1)
%     baseline_avg(i,j)=mean(squeeze(F_tan(i,end-700:end-100,j)));
%     std_dev(i,j)=std(squeeze(F_tan(i,end-700:end-100,j)));
% end
% j=2;
% for i=1:(cou_nt-1)
%     baseline_avg(i,j)=mean(squeeze(F_tan(i,100:700,j)));
%     std_dev(i,j)=std(squeeze(F_tan(i,100:700,j)));
% end

figure
j=1
for i=1:cou_nt-1
    plot(def_ang(i,:,j),F_tan(i,:,j),'-','Color',col(i,:))
    hold on
end
yline(0)
j=2
for i=1:cou_nt-1
    plot(def_ang(i,:,j),F_tan(i,:,j),'--','Color',col(i,:))
    hold on
end
% yline(0)
hold off
xlabel("$\theta$ [mrad]",'Interpreter','latex')
ylabel("$F_v$ [N]",'Interpreter','latex')



%% Jose's modified backpropagation for signal isolation

for j=1:2
    for i=1:cou_nt-1
        sig=squeeze(F_tan(i,:,j)); %N
        sig_x=squeeze(def_ang(i,:,j)); %mrad
        ref_x=baseline_avg(i,j); %N
        [max_torque_org(i,j),max_torque_pos_org(i,j)]=max(F_tan(i,:,j));
        [min_torque_org(i,j),min_torque_pos_org(i,j)]=min(F_tan(i,:,j));
        sig=abs(sig);
        backprop_window=6; %number of points to use for slope calculation
        slope_offset=100; % offeset to the left of max or min position
        if j==1
            for k=2:max_torque_pos_org(i,j)-slope_offset-backprop_window-1
            end_pos=max_torque_pos_org(i,j)-slope_offset-k;
            start_pos=end_pos-backprop_window;
            fit_lin=polyfit(sig_x(start_pos:end_pos),sig(start_pos:end_pos),1);
            m(k)=fit_lin(1);
            start_array(k)=start_pos;
            end_array(k)=end_pos;
            y_int(k)=fit_lin(2);
            x_intersect=(ref_x-fit_lin(2))/fit_lin(1);
            x_int(k)=x_intersect;
            x_start(k)=sig_x(start_pos);
            x_end(k)=sig_x(end_pos);
            if and(x_int(k)>x_end(k),x_int(k-1)<x_start(k-1))
                lim_start_index(i,j)=start_array(k);
                break
            end
            end
    
%           figure
%           plot(x_int)
%           hold on
%           plot(x_start)
%           plot(x_end)    
%           hold off
    
    clear end_pos start_pos fit_lin m start_array end_array y_int x_intersect x_int x_start x_end k
    
  
    sig=fliplr(sig);
    sig_x=fliplr(sig_x);
    s=size(sig);
    peakpos=s(2)-min_torque_pos_org(i,j);
    for k=2:peakpos-slope_offset-backprop_window-1
        end_pos=peakpos-slope_offset-k;
        start_pos=end_pos-backprop_window;
        fit_lin=polyfit(sig_x(start_pos:end_pos),sig(start_pos:end_pos),1);
        m(k)=fit_lin(1);
        start_array(k)=start_pos;
        end_array(k)=end_pos;
        y_int(k)=fit_lin(2);
        x_intersect=(ref_x-fit_lin(2))/fit_lin(1);
        x_int(k)=x_intersect;
        x_start(k)=sig_x(start_pos);
        x_end(k)=sig_x(end_pos);
        if and(x_int(k-1)>x_start(k-1),x_int(k)<x_end(k))
            lim_stop_index(i,j)=s(2)-end_array(k);
            break
        end
        
    end
    
% figure
%     plot(x_int)
%     hold on
%     plot(x_start)
%     plot(x_end)   
%     hold off
    
    clear end_pos start_pos fit_lin m start_array end_array y_int x_intersect x_int x_start x_end k
else
    
   for k=2:min_torque_pos_org(i,j)-slope_offset-backprop_window-1
        end_pos=min_torque_pos_org(i,j)-slope_offset-k;
        start_pos=end_pos-backprop_window;
        fit_lin=polyfit(sig_x(start_pos:end_pos),sig(start_pos:end_pos),1);
        m(k)=fit_lin(1);
        start_array(k)=start_pos;
        end_array(k)=end_pos;
        y_int(k)=fit_lin(2);
        x_intersect=(ref_x-fit_lin(2))/fit_lin(1);
        x_int(k)=x_intersect;
        x_start(k)=sig_x(start_pos);
        x_end(k)=sig_x(end_pos);
        if and(x_int(k)<x_end(k),x_int(k-1)>x_start(k-1))
            lim_start_index(i,j)=end_array(k);
            break
        end
    end
    
%     figure
%     plot(x_int)
%     hold on
%     plot(x_start)
%     plot(x_end)    
%     hold off
    
    clear end_pos start_pos fit_lin m start_array end_array y_int x_intersect x_int x_start x_end k
    
    
    sig=fliplr(sig);
    sig_x=fliplr(sig_x);
    s=size(sig);
    peakpos=s(2)-max_torque_pos_org(i,j);
    for k=2:peakpos-slope_offset-backprop_window-1
        end_pos=peakpos-slope_offset-k;
        start_pos=end_pos-backprop_window;
        fit_lin=polyfit(sig_x(start_pos:end_pos),sig(start_pos:end_pos),1);
        m(k)=fit_lin(1);
        start_array(k)=start_pos;
        end_array(k)=end_pos;
        y_int(k)=fit_lin(2);
        x_intersect=(ref_x-fit_lin(2))/fit_lin(1);
        x_int(k)=x_intersect;
        x_start(k)=sig_x(start_pos);
        x_end(k)=sig_x(end_pos);
        if and(x_int(k-1)<x_start(k-1),x_int(k)>x_end(k))
            lim_stop_index(i,j)=s(2)-start_array(k);
            break
        end
        
    end
     
% figure
%     plot(x_int)
%     hold on
%     plot(x_start)
%     plot(x_end)    
%     hold off
    
    clear end_pos start_pos fit_lin m start_array end_array y_int x_intersect x_int x_start x_end k
    


end
    end
end


for j=1:2
    for i=1:cou_nt-1
%         figure        
%         plot(def_ang(i,:,j),F_tan(i,:,j))
%         hold on
%         plot(def_ang(i,lim_start_index(i,j):lim_stop_index(i,j),j),F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j))
            angular_displacement(i,j)=(def_ang(i,lim_start_index(i,j),j)-def_ang(i,lim_stop_index(i,j),j))/1000; 
%         hold off
%         xlabel("Deflection angle (mrad)")
%         ylabel("Tangential Forces (N)")
    end
    
end
%% check_limits
for i=1:cou_nt-1
    for j=1:2
        def_angle_start(i,j)=def_ang(i,lim_start_index(i,j),j);
        def_angle_end(i,j)=def_ang(i,lim_stop_index(i,j),j);
    end
end
contact_def_angle_geometric=mean(def_angle_start(:,2))-mean(def_angle_start(:,1));
mean_contact_start_def_angle=[mean(def_angle_start(:,1)) mean(def_angle_start(:,2))];
mean_contact_mid_def_angle=[mean_contact_start_def_angle(1)+(contact_def_angle_geometric/2) mean_contact_start_def_angle(2)-(contact_def_angle_geometric/2)];
mean_mid_def_angle=mean(mean_contact_mid_def_angle);
dist_from_mid_defangle=abs(def_ang-mean_mid_def_angle);
[val_min,val_loc]=min(dist_from_mid_defangle(:,:,1),[],2);
[val_min(:,2),val_loc(:,2)]=min(dist_from_mid_defangle(:,:,2),[],2);
contact_def_len_geometric=contact_def_angle_geometric*moment_arm/1000;
contact_def_len_geometric_mean=mean(contact_def_len_geometric);
contact_def_angle_len(1:cou_nt-1,1)=def_angle_end(:,1)-mean(def_angle_start(:,1));
contact_def_angle_len(1:cou_nt-1,2)=mean(def_angle_start(:,2))-def_angle_end(:,2);
contact_def_len=contact_def_angle_len*moment_arm/1000;
for j=1:2
    for i=1:cou_nt-1
        ran_ge=val_loc(i,j)-5:val_loc(i,j)+5;
        ran_ge_vert=val_loc(i,j)-20:val_loc(i,j)+20;
        fit_defang_mid=polyfit(def_ang(i,ran_ge,j),F_tan(i,ran_ge,j),1);
        F_tan_mid(i,j)=polyval(fit_defang_mid,mean_mid_def_angle);
        fit_vert_force_mid=polyfit(def_ang(i,ran_ge_vert,j),F_norm(i,ran_ge_vert,j),2);
        F_vert_mid(i,j)=polyval(fit_vert_force_mid,mean_mid_def_angle);
        v_hz_mean(i,j)=mean(v_hz(i,:,j));
        mu_exp(i,j)=F_tan_mid(i,j)/F_vert_mid(i,j);
        
%         figure
%         scatter(def_ang(i,ran_ge,j),F_tan(i,ran_ge,j))
%         hold on
%         plot(def_ang(i,ran_ge,j),polyval(fit_defang_mid,def_ang(i,ran_ge,j)))
%         hold off
%         xlabel('defangle')
%         ylabel('tangential force')
%         title('central force tan fit')
%         
%         figure
%         scatter(def_ang(i,ran_ge_vert,j),F_norm(i,ran_ge_vert,j))
%         hold on
%         plot(def_ang(i,ran_ge_vert,j),polyval(fit_vert_force_mid,def_ang(i,ran_ge_vert,j)))
%         hold off
%         xlabel('defangle')
%         ylabel('tangential force')
%         title('central force norm fit')
    end
end


%% forces at center of contact
figure
plot(v_hz_mean(:,1),F_tan_mid(:,1))
hold on
plot(abs(v_hz_mean(:,2)),abs(F_tan_mid(:,2)))
hold off
xlabel('v [m/s]')
ylabel('F_h^{mid}')

figure
plot(v_hz_mean(:,1),F_vert_mid(:,1))
hold on
plot(abs(v_hz_mean(:,2)),abs(F_vert_mid(:,2)))
hold off
xlabel('v [m/s]')
ylabel('F_v^{mid}')

figure
plot(v_hz_mean(:,1),mu_exp(:,1))
hold on
plot(abs(v_hz_mean(:,2)),abs(mu_exp(:,2)))
hold off
xlabel('v [m/s]')
ylabel('\mu^{mid}')



%% contact length vs. velocity
figure
plot(v_hz_mean(:,1),contact_def_len(:,1))
hold on
plot(abs(v_hz_mean(:,2)),contact_def_len(:,2))
hold off
xlabel('v [m/s]')
ylabel('l_c [m]')

%% mean of response parameters
dim_array=size(v_hz_mean)
for i=1:dim_array(1)
    F_tan_mid_avg(i)=mean(F_tan_mid(i,:));
    F_vert_mid_avg(i)=mean(F_vert_mid(i,:));
    mu_exp_mean(i)=mean(abs(mu_exp(i,:)));
    v_hz_mean_avg(i)=mean(abs(v_hz_mean(i,:)));
    contact_def_len_mean(i)=mean(abs(contact_def_len(i,:)));
end

%% elastic displacement
for i=1:cou_nt-1
    max_ang_disp(i)=(def_ang(i,lim_start_index(i,1),1)-def_ang(i,lim_start_index(i,2),2))/1000; %mrad converted to rad
end
max_lin_disp_m=abs(max_ang_disp*moment_arm);
max_lin_disp_m_avg= mean(max_lin_disp_m);

max_lin_disp_mm=abs(max_ang_disp*moment_arm*1000);
lin_displacement=abs(angular_displacement)*moment_arm;
lin_displacement_mm=lin_displacement*1000;

%% seperation between point of contact and detachment
for i=1:cou_nt-1
    processed_data.sep_clock(i,1)=moment_arm*(def_ang(i,lim_start_index(i,2),2)-def_ang(i,lim_stop_index(i,1),1))/1000;
    processed_data.sep_clock(i,2)=moment_arm*(def_ang(i,lim_start_index(i,1),1)-def_ang(i,lim_stop_index(i,2),2))/1000;
    processed_data.max_amp(i,1)=max(F_tan(i,:,1));
    processed_data.min_amp(i,1)=min(F_tan(i,:,1));
    processed_data.max_amp(i,2)=max(-F_tan(i,:,2));
    processed_data.min_amp(i,2)=min(-F_tan(i,:,2));
end

%% contact length in one direction
processed_data.contact_len=moment_arm*angular_displacement/1000;
processed_data.avg_contact_len=(abs(processed_data.contact_len(:,1))+abs(processed_data.contact_len(:,2)))/2;

%% point crossing in time
figure
for i=1:cou_nt-1
j=1;
plot(time_exp(i,lim_start_index(i,j):lim_stop_index(i,j),j)-time_exp(i,lim_start_index(i,j),j),F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg(i,j),"Color",col(i,:),"LineStyle","-")
hold on
j=2;
plot(time_exp(i,lim_start_index(i,j):lim_stop_index(i,j),j)-time_exp(i,lim_start_index(i,j),j),F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg(i,j),"Color",col(i,:),"LineStyle","-.")
yline(0)
end
hold off
xlabel("time (s)")
ylabel("Torque(mN.m)")

%% point crossing in absolute displacement

%color matrix
col_mat_green=[152, 241, 217;
    125, 195, 156;
    103, 149, 100;
    82, 105, 51;
    59, 64, 4]/255;
col_mat_purple=[161, 113, 184;
    135, 76, 135;
    105, 41, 89;
    73, 3, 48]/255;

figure
for i=1:cou_nt-1
j=1;
name_tag=strcat(num2str(mean(v_hz(i,:,j))), " m/s");

plot(abs(def_ang(i,lim_start_index(i,j):lim_stop_index(i,j),j)-def_ang(i,lim_start_index(i,j),j))*moment_arm/1000,F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg(i,j),"Color",col(i,:),"LineStyle","-","DisplayName",name_tag)

% plot(abs(def_ang(i,lim_start_index(i,j):lim_stop_index(i,j),j))*moment_arm/1000,F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg(i,j),"Color",col_mat_purple(i,:),"LineStyle","-")
processed_data.tofit(i,j).x=(squeeze(abs(def_ang(i,lim_start_index(i,j):lim_stop_index(i,j),j)-def_ang(i,lim_start_index(i,j),j))))*(moment_arm/1000); % divide by 1000 to change mrad to rad
processed_data.tofit(i,j).F_hz=squeeze(F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg(i,j));
processed_data.tofit(i,j).ver=squeeze(F_norm(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg_vert(i,j));
processed_data.tofit(i,j).hz_by_ver=processed_data.tofit(i,j).F_hz./processed_data.tofit(i,j).ver;
processed_data.tofit(i,j).time_array=time_exp(i,lim_start_index(i,j):lim_stop_index(i,j),j)-time_exp(i,lim_start_index(i,j),j);
processed_data.tofit(i,j).v_overlap=v_hz(i,lim_start_index(i,j):lim_stop_index(i,j),j); % velocity of particle particle shear

hold on
j=2;

% plot(abs(def_ang(i,lim_start_index(i,j):lim_stop_index(i,j),j)-def_ang(i,lim_start_index(i,j),j)),F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg(i,j),"Color",col(i,:),"LineStyle","-")


% plot(abs(def_ang(i,lim_start_index(i,j):lim_stop_index(i,j),j))*moment_arm/1000,F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg(i,j),"Color",col_mat_green(i,:),"LineStyle","-")
processed_data.tofit(i,j).x=(squeeze(abs(def_ang(i,lim_start_index(i,j):lim_stop_index(i,j),j)-def_ang(i,lim_start_index(i,j),j))))*(moment_arm/1000); % divide by 1000 to change mrad to rad
processed_data.tofit(i,j).F_hz=squeeze(F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg(i,j));
processed_data.tofit(i,j).ver=squeeze(F_norm(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg_vert(i,j));
processed_data.tofit(i,j).hz_by_ver=processed_data.tofit(i,j).F_hz./processed_data.tofit(i,j).ver;
processed_data.tofit(i,j).hz_by_ver=processed_data.tofit(i,j).F_hz./processed_data.tofit(i,j).ver;
processed_data.tofit(i,j).hz_by_ver=processed_data.tofit(i,j).F_hz./processed_data.tofit(i,j).ver;

processed_data.tofit(i,j).time_array=time_exp(i,lim_start_index(i,j):lim_stop_index(i,j),j)-time_exp(i,lim_start_index(i,j),j);
processed_data.tofit(i,j).v_overlap=v_hz(i,lim_start_index(i,j):lim_stop_index(i,j),j); % velocity of particle particle shear



end
yline(0,"DisplayName","Baseline")
hold off
xlabel("Displacement (m)")
ylabel("Horizontal Force(N)")
lgd = legend('Location','northeast');
lgd.Title.String = 'Rotational Speed';

figure
for i=1:cou_nt-1
j=1;

plot((abs(def_ang(i,lim_start_index(i,j):lim_stop_index(i,j),j)-def_ang(i,lim_start_index(i,j),j))*moment_arm/1000)/(mean([particle_dia_lower_perp particle_dia_upper_perp])/2),(F_tan(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg(i,j))/overlap_avg(i,j),"Color",col(i,:),"LineStyle","-")
hold on
end
yline(0)
hold off
xlabel('^{Displacement}/_{Avg. Radius}')
ylabel('^{Horizontal Force(N)}/_{Max. overlap(m)}')



   %% plotting vertical forces at different levels 
figure
for i=1:cou_nt-1
j=1;
name_tag=strcat(num2str(mean(v_hz(i,:,j))), " mrad/s");
plot(processed_data.tofit(i,j).x,processed_data.tofit(i,j).ver,"Color",col(i,:),"LineStyle","-","DisplayName",name_tag);
hold on
% j=2;
% plot(processed_data.tofit(i,j).x,processed_data.tofit(i,j).ver,"Color",col(i,:),"LineStyle","-.");

end
yline(0,"DisplayName","Baseline")
hold off
xlabel("Absolute Angular displacement (m)")
ylabel("Vertical Forces (N)")
lgd = legend('Location','south');
lgd.Title.String = 'Rotational Speed';

%% zero crossing

for j=1:2
for i=1:cou_nt-1
ref_x=baseline_avg(i,j);
[max_torque_org(i,j),max_torque_pos_org(i,j)]=max(F_tan(i,:,j));
[min_torque_org(i,j),min_torque_pos_org(i,j)]=min(F_tan(i,:,j));
window_def=45*std_dev(i,j);
vec_range=min(min_torque_pos_org(i,j),max_torque_pos_org(i,j)):max(min_torque_pos_org(i,j),max_torque_pos_org(i,j));
%left_lim
logic_up=F_tan(i,vec_range,j)<ref_x+window_def;
logic_down=F_tan(i,vec_range,j)>ref_x-window_def;
if j==1
    start_pos=find(logic_up,1,'first');
    end_pos=find(logic_down,1,'last');
else
    start_pos=find(logic_down,1,'first');
    end_pos=find(logic_up,1,'last');
end

% %plot to check limits
% figure
% plot(F_tan(i,vec_range,j),"DisplayName","data")
% hold on
% yline(ref_x,"DisplayName","reference")
% yline((ref_x+window_def),"DisplayName","upper limit","Color","Red")
% yline((ref_x-window_def),"DisplayName","lower limit","Color","Green")
% xline(start_pos,"Color","Red","LineStyle","--")
% xline(end_pos,"Color","Green","LineStyle","--")
% hold off
% legend
% ylabel('Hz Forces')
% xlabel('indices')
% title("check limits")

cross_fit=polyfit(def_ang(i,vec_range(1)+(start_pos:end_pos)-1,j),F_tan(i,vec_range(1)+(start_pos:end_pos)-1,j),1);
cross_pt(i,j)=(ref_x-cross_fit(2))/cross_fit(1);
end
end

% for j=1:2
% for i=1:4
% 
% figure
% plot(def_ang(i,:,j),F_tan(i,:,j))
% hold on
% plot(cross_pt(i,j),ref_x,'r*')
% yline(ref_x)
% hold off
% xlabel("Deflection angle (mrad)")
% ylabel("Torque (mN.m)")
% end
% end


%% fit normal force data

% skewed gaussian curve
i=1;
j=1;
x=squeeze(def_ang(i,:,j));
y=squeeze(F_norm(i,:,j)- baseline_avg_vert(i,j));
sd=std(x*6);
mod.gauss=skewed_gaussian(x,y,sd);

% 2nd degree polynomial
i=1;
j=1;
x=squeeze(def_ang(i,lim_start_index(i,j):lim_stop_index(i,j),j));
y=squeeze(F_norm(i,lim_start_index(i,j):lim_stop_index(i,j),j)- baseline_avg_vert(i,j));
fit_poly=polyfit(x,y,2);
% figure
% scatter(x,y)
% hold on
% plot(x,polyval(fit_poly,x))
% hold off


% figure
% plot(x,polyval(fit_poly,x)-y)



%% histogram of averaged rate of torque change distribution data

window_hist=300;

%clockwise
%direction levels
% figure
%     for i=1:(cou_nt-1)
%         subplot(2,2,i)
%         hist(squeeze(dtorque_dtheta_movmean(i,:,1)),window_hist)
%         xlabel("dtorque/dtheta movmean")
%         ylabel("frequency")
%     end    
% sgtitle("dTorque/dDeflection angle Rate Distribution Clockwise")


for i=1:(cou_nt-1)
    [bin(i).N(1,:),bin(i).edges(1,:)] = histcounts(squeeze(dtorque_dtheta_movmean(i,:,1)),window_hist);
    [mode,I]=max(bin(i).N(1,:));
    int_start_dTorddef(1,i)=bin(i).edges(1,I);
    int_prev_dTorddef(1,i)=bin(i).edges(1,I-1);
    int_end_dTorddef(1,i)=bin(i).edges(1,I+1);
end

% anticlockwise
% figure
% for i=1:(cou_nt-1)
%     subplot(2,2,i)
%     hist(squeeze(dtorque_dtheta_movmean(i,:,2)),window_hist)
%     xlabel("dtorque/dtheta movmean")
%     ylabel("frequency")
% end
% sgtitle('dTorque/dDeflection angle Rate Distribution Anticlockwise')

for i=1:(cou_nt-1)
    [bin(i).N(2,:),bin(i).edges(2,:)] = histcounts(squeeze(dtorque_dtheta_movmean(i,:,2)),window_hist);
    [mode,I]=max(bin(i).N(2,:));
    int_start_dTorddef(2,i)=bin(i).edges(2,I);
    int_prev_dTorddef(2,i)=bin(i).edges(2,I-1);
    int_end_dTorddef(2,i)=bin(i).edges(2,I+1);
end

%% torque isolatation by dtorque/dtheta 
len_veccheck=100;

for j=1:2 % direction of rotation
%     figure
    for i=1:(cou_nt-1) % speed level
%  
%     subplot(2,2,i)
    [max_torque(i,j),max_torque_pos(i,j)]=max(F_tan(i,:,j));
    [min_torque(i,j),min_torque_pos(i,j)]=min(F_tan(i,:,j));
    size_torque=size(squeeze(F_tan(i,:,j)));
    logic_vec=ones(size_torque);
    logic_vec(1)=0;
    logic_vec(1:max_torque_pos(i,j))=squeeze(~and(dtorque_dtheta_movmean(i,1:max_torque_pos(i,j),j)<int_end_dTorddef(j,i),dtorque_dtheta_movmean(i,1:max_torque_pos(i,j),j)>int_start_dTorddef(j,i)));
    logic_vec(min_torque_pos(i,j)+1:end)=squeeze(~and(dtorque_dtheta_movmean(i,min_torque_pos(i,j):end,j)<int_end_dTorddef(j,i),dtorque_dtheta_movmean(i,min_torque_pos(i,j):end,j)>int_start_dTorddef(j,i)));
    logic_vec_right=zeros(size(logic_vec));
    logic_vec_left=zeros(size(logic_vec));
    for k=1:length(logic_vec)-len_veccheck
        if logic_vec(k:k+len_veccheck)==1;
            logic_vec_right(k)=1;
            index_start(i,j)=k;
            break;
        end
    end
    for k=length(logic_vec):-1:1+len_veccheck
        if logic_vec(k-len_veccheck:k)==1;
            logic_vec_left(k)=1;
            index_end(i,j)=k;
            break;
        end
    end
%     plot(logic_vec)
%     hold on
%     plot(logic_vec_right)
%     hold on
%     plot(logic_vec_left)
%     
%     plot(F_tan(i,:,j)./max(F_tan(i,:,j)))
    
%     plot(dtorque_dtheta_movmean(i,:,j)./max(dtorque_dtheta_movmean(i,:,j)))
%     
%     plot(dtorque_dtheta(i,:,j)./max(dtorque_dtheta(i,:,j)))
%     hold off
    
    
    
    end
end
% sgtitle('signal isolation by torque histogram')

%% checking clipped data
for i=1:cou_nt-1
%     figure
        x1=def_ang(i,index_start(i,1):index_end(i,1),1);
        y1=F_tan(i,index_start(i,1):index_end(i,1),1)-baseline_avg(i,1);
        x2=def_ang(i,index_start(i,2):index_end(i,2),2);
        y2=F_tan(i,index_start(i,2):index_end(i,2),2)-baseline_avg(i,2);
        min_defang=min([x1 x2]);
        max_defang=max([x1 x2]);
        x1=(x1-min_defang)*moment_arm;
        x2=(x2-min_defang)*moment_arm;
        total_disp=(max_defang-min_defang)*moment_arm;
%         plot(x1,y1)
%         hold on
%         plot(x2,y2)
%         yline(0)
%         hold off
end


%% velocity averaging
for j=1:2
    for i=1:cou_nt-1
    processed_data.tofit(i,j).v_hz_mean=abs(mean(v_hz(i,lim_start_index(i,j):lim_stop_index(i,j),j)));
    processed_data.v_hz_mean(i,j)=processed_data.tofit(i,j).v_hz_mean;
    end

end

% processed_data.tofit(i,j).v_hz_mean=v_hz_mean;

%% passing experimental crossing point
for j=1:2
    for i=1:cou_nt-1
        processed_data.tofit(i,j).x_zero_exp=abs((cross_pt(i,j)-def_ang(i,lim_start_index(i,j),j))/1000)*moment_arm;
        processed_data.x_zero_exp(i,j)=processed_data.tofit(i,j).x_zero_exp;
    end
end

% processed_data.tofit(i,j).x_zero_exp=x_zero_exp;
 processed_data.tofit(i,j).overlap_avg(1:cou_nt-1,1:2)=overlap_avg;
 
%% calculating overlap numerically
processed_data.overlap_max_numeric=overlap_max_r1_r2_true(particle_dia_lower_perp/2,particle_dia_upper_perp/2,max_lin_disp_m_avg,moment_arm);



%% returning values to main function
processed_data.F_hz=F_tan;
processed_data.F_ver=F_norm;
processed_data.F_ver_o=F_norm_o;
processed_data.F_norm_max=F_norm_max;
processed_data.F_hz_mid=F_tan_mid;
processed_data.F_vert_mid=F_vert_mid;
processed_data.contact_def_len=contact_def_len;
processed_data.contact_def_len_geometric=contact_def_len_geometric;
processed_data.v_hz_mean=v_hz_mean;
processed_data.time=time_exp;
processed_data.def_ang=def_ang;
processed_data.overlap=overlap_avg;
processed_data.F_tan_mid_avg=F_tan_mid_avg;
processed_data.F_vert_mid_avg=F_vert_mid_avg;
processed_data.mu_exp_mean=mu_exp_mean;
processed_data.v_hz_mean_avg=v_hz_mean_avg;
processed_data.contact_def_len_mean=contact_def_len_mean;
processed_data.def_angle_start=def_angle_start;
% processed_data.v_hz_mean=v_hz_mean;
processed_data.lin_disp_m=max_lin_disp_m;
processed_data.stat.baseline_avg_hz=baseline_avg;
processed_data.stat.baseline_avg_vert=baseline_avg_vert;
processed_data.stat.std_dev_hz=std_dev;
processed_data.stat.baseline_avg_vert=baseline_avg_vert;
processed_data.filter.lim_start_index=lim_start_index;
processed_data.filter.lim_stop_index=lim_stop_index;
processed_data.filter.angular_displacement=angular_displacement;
processed_data.filter.lin_displacement=lin_displacement;
processed_data.filter.cross_pt_defangle=cross_pt;
processed_data.filter.max_lin_disp_m_avg=max_lin_disp_m_avg;

end
