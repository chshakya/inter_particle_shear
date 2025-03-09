clear all
close all

cd '/home/ubuntu_pcc3/Downloads/matlab_script-examples'

Ge=0.4e4;
G1=5*Ge;
vp=0.5;
T=2000;
D_real=mean([18.16,18.55,18.79,19.32])
R=(D_real/2)/1000;
str_ain=0.018;
% dmax=str_ain*2*2*R;
dmax=6.8e-4
mu=5e-3;

data_set(1).data=viscoelastic_experiment(Ge,G1,vp,T,R,dmax(1),mu);

for i=2:5
dmax(i)=dmax(i-1)+(0.1/1000);
data_set(i).data=viscoelastic_experiment(Ge,G1,vp,T,R,dmax(i),mu);
end

y_theo=2*R-dmax;
contact_length=(((2*R)^2)-y_theo.^2).^0.5;
offset_x=(contact_length(end)-contact_length);

f=figure
for i=1:5
    name_tag=sprintf('%0.2e',dmax(i))
    scatter(data_set(i).data.v,data_set(i).data.x_zero./data_set(i).data.x_geometric,'DisplayName',name_tag)
    hold on
end
hold off
xlabel("$v$ [m/s]",'Interpreter','latex')
ylabel("$\frac{l_c}{l_g}$",'Interpreter','latex')
lgd = legend('Location','southeast','Interpreter','latex');
lgd_title='$\delta_{max} [m]$';
lgd.Title.String = lgd_title;
lgd.FontSize = 18;
ax1 = gca;
ax1.LineWidth=2;
ax1.FontSize = 20;
ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')
f.Position=[200 300 560 500]
AddLetters2Plots({ax1},{'b)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,'/home/ubuntu_pcc3/paper_writing/inter_particle_shear/l_c_by_lg_sim.svg')
saveas(gcf,'/home/ubuntu_pcc3/paper_writing/inter_particle_shear/l_c_by_lg_sim','epsc')

figure
for i=1:4
    plot(data_set(1).data.x_plot(i,:),-data_set(1).data.Resx(i,:))
    hold on
end
yline(0)
xlabel('Displacement[m]')
ylabel('F_{hz} [N]')
hold off
saveas(gcf,'/home/ubuntu_pcc3/paper_writing/inter_particle_shear/F_hz_ver_disp_sim.svg')
saveas(gcf,'/home/ubuntu_pcc3/paper_writing/inter_particle_shear/F_hz_ver_disp_sim','epsc')
saveas(gcf,'/home/ubuntu_pcc3/paper_writing/inter_particle_shear/F_hz_ver_disp_sim.pdf')

f=figure
for i=1:5
    name_tag=sprintf('%0.2e',dmax(i))
    plot(4.95e-3+data_set(i).data.x_plot(1,:)+offset_x(i),-data_set(i).data.Resx(1,:),'DisplayName',name_tag,'LineWidth',2)
    hold on
end
yline(0,'HandleVisibility','off')
hold off
xlim([4.98e-3 17e-3])
xlabel("$\theta$ [m/s]",'Interpreter','latex')
ylabel("$F_h$ [N]",'Interpreter','latex')
lgd = legend('Location','southwest','Interpreter','latex');
lgd_title='$\delta_{max} [m]$';
lgd.Title.String = lgd_title;
lgd.FontSize = 16;
ax1 = gca;
ax1.LineWidth=2;
ax1.FontSize = 20;
ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')
f.Position=[200 300 560 500]
AddLetters2Plots({ax1},{'c)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,'/home/ubuntu_pcc3/paper_writing/inter_particle_shear/hz_force_diff_overlap_sim.svg')
saveas(gcf,'/home/ubuntu_pcc3/paper_writing/inter_particle_shear/hz_force_diff_overlap_sim','epsc')


f=figure
for i=1:5
    name_tag=sprintf('%0.2e',dmax(i))
    scatter(data_set(i).data.v,data_set(i).data.x_zero./data_set(i).data.x_geometric,'DisplayName',name_tag,'LineWidth',3)
    hold on
end
hold off
ylim([0.94 1])
box on
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
f.Position=[200 300 560 500]
AddLetters2Plots({ax1},{'d)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,'/home/ubuntu_pcc3/paper_writing/inter_particle_shear/l_c_by_lg_sim.svg')
saveas(gcf,'/home/ubuntu_pcc3/paper_writing/inter_particle_shear/l_c_by_lg_sim','epsc')


function analytical=viscoelastic_experiment(Ge,G1,vp,T,R,dmax,mu)
y=2*R-dmax;  %vertical distance between centers of spheres
Rstar=R/2;

v_char=sqrt(4*R^2-y^2)/T;
analytical.v_char=v_char;


v_fact=2.38e-5*[1 2 5 10];%221);
col=hsv(2*length(v_fact));




for m=1:length(v_fact)
% m=1;
v=v_fact(m);%*3.96e-04;  %velocity of moving sphere
 
%define t=to and x=0 as point where spheres are on top of each other and 
%time of first contact:d=2*R
 
to=sqrt(4*R^2-y^2)/v;  %time at which beads are on top of each other
dt=to/200;  %time increment
t=[0:dt:2*to];   
x=v*(t); %x-position of sphere
xo=v*to;%x-position of fixed sphere
d=sqrt(y^2+(x-xo).^2); 
del=(2*R-d); %indentation
%  
% ddeldt=-v^2*(t-to)./(2*d); %ddelta/dt=-0.5dd/dt=-(x/2d)*v
 
%Elastic Hertzian part:
Fe=(4*Ge*sqrt(Rstar)/(3*(1-vp)))*del.^(3/2);  %elastic force

analytical.force_struct(m).x=x;
analytical.force_struct(m).Fe=Fe;
analytical.force_struct(m).del=del;
analytical.force_struct(m).Fe_hz=(x-xo).*Fe./d;

%dissipation during elastic contact (ideally=0) calculated here to look at
%baseline numeric noise
dissipation_elastic_x(m)=trapz(x,((x-xo).*Fe./d));
 
%Dissipative part:
pref=4*sqrt(Rstar)/(2*(1-vp));
ds=dt/200;  %time step for numerical integration (eq 23)
Fd=zeros(size(t));
 
for i=2:length(t)
    s=[0:ds:t(i)];
    xi=v*(s);
    di=sqrt(y^2+(xi-xo).^2);
    deli=(2*R-di);
    d12i=deli.^(1/2);
    ddoti=-(xi-xo)*v./(2*di);
    del_t=t(i)-s;
    expo=exp(-(del_t)/T);
    integ(i).arr=G1*exp(-(del_t)/T).*d12i.*ddoti;
    Fd(i)=pref*trapz(s,integ(i).arr); %numerical integral for Fd
end

analytical.force_struct(m).Fd=Fd;
% Fd=Fd.*(Fd>0);

% figure
% plot(Fd)
% hold on
% yline(0)
% hold off

Res(m,:)=Fe+Fd;
analytical.force_struct(m).Fe_Fd_raw=Fe+Fd;

trigger_val=0; %set trigger value to 0 assuming there ae no tensile forces
x_zero(m)=x(end);
%find cycle endpoint between displacement steps(x-field)
[max_Res,max_Res_pos]=max(Res(m,:));
[min_Res,min_Res_pos]=min(Res(m,:));
for i=min(min_Res_pos,max_Res_pos)+10:max(min_Res_pos,max_Res_pos)-1
    if or(and(Res(m,i)<0,Res(m,i+1)>0),and(Res(m,i)>0,Res(m,i+1)<0))
        zero_fit=polyfit(x(i:i+1),Res(m,i:i+1),1);
        x_zero_array=x(i:i+1);
        y_zero_array=Res(m,i:i+1);
        x_zero(m)=-zero_fit(2)/zero_fit(1);
        trigger_val=1; %set trigger value to 1 to kill tensile forces
        break
    end
end

%formation of new arrays based on analytical cycle end (interpolated end x)
if trigger_val==1 % killing tensile forces
    x_cut(m).contact_length=[(x(1:i)),(x_zero(m))];
    res_cut(m).resultant_viscoelastic=[Res(m,1:i),0];
    res_cut(m).resultant_viscoelastic_x=[(x(1:i)-xo).*Res(m,1:i)./d(1:i),0];
else
    x_cut(m).contact_length=x;
    res_cut(m).resultant_viscoelastic=Res(m,:);
    res_cut(m).resultant_viscoelastic_x=(x-xo).*Res(m,:)./d;
end

dissipation_viscous_hertz(m)=trapz((x_cut(m).contact_length),(res_cut(m).resultant_viscoelastic_x));

analytical.force_struct(m).dissipation_viscous_hertz=dissipation_viscous_hertz(m);
analytical.force_struct(m).resultant_visoelastic_eff=res_cut(m).resultant_viscoelastic;
analytical.force_struct(m).resultant_visoelastic_eff_x=res_cut(m).resultant_viscoelastic_x;

%% friction dissipation

res_cut(m).friction=mu*res_cut(m).resultant_viscoelastic;
if trigger_val==1
res_cut(m).friction_x= [mu*y*Res(m,1:i)./d(1:i),0];%y*res_cut(m).friction(m,:)./[d(1:i),0];
else
res_cut(m).friction_x= mu*y*Res(m,:)./d;%y*res_cut(m).friction(m,:)./[d(1:i),0];   
end

dissipation_friction(m)=trapz(x_cut(m).contact_length,res_cut(m).friction_x);
analytical.force_struct(m).dissipation_friction=dissipation_friction(m);

%formation of new resultant force array based on discrete x
for i=1:length(Res(m,:))
    if and(Fe(i)>0,Res(m,i)>0)
        Res_mod(m,i)=Res(m,i);
    elseif and(Fe(i)<0,Res(m,i)<0)
        Res_mod(m,i)=Res(m,i);
    else
       Res_mod(m,i)=0;   
    end
end

Array=Res_mod(m,2:end);
cross_end_index(m)=find(Array==0,1,'first')+1;
x_end(m)=x(cross_end_index(m));
dist_char(m)=v*T;

%Frictional Part
FF=mu*Res_mod(m,:);


FxH(m,:)=(x-xo).*Fe./d;  %x-component for force, Hertz contribution
FyH(m,:)=y.*Fe./d;  %y-component for force, Hertz contribution
FxD(m,:)=(x-xo).*Fd./d;  %x-component for force, dissipative contribution
FyD(m,:)=y.*Fd./d;  %y-component for force, dissipative contribution

Resx(m,:)= (x-xo).*Res_mod(m,:)./d;
Resy(m,:)= y.*Res_mod(m,:)./d;
FxF(m,:)=y.*FF./d;  %x-component for force, Frictional contribution
FyF(m,:)=(x-xo).*FF./d; %y-component for force, Frictional contribution



x_plot(m,:)=x;
x_geometric=x(end);

end

%% contact length vs. velocity
% figure
% scatter(v_fact,x_zero)

%% minimum cycle length
normalized_x=x_zero/x(end);
normalized_vel=v_fact/v_char;
[x_min,x_min_pos]=min(normalized_x);
vel_min_x_norm=normalized_vel(x_min_pos);


%% data to pass
analytical.x_norm=normalized_x;
analytical.v_norm=normalized_vel;
analytical.x_min=x_min;
analytical.vel_at_min_x_norm=vel_min_x_norm;
analytical.v=v_fact;
analytical.x_zero=x_zero;
analytical.x_geometric=x_geometric;
analytical.v_char=v_char;
analytical.x_plot=x_plot;
analytical.Resx=Resx-FxF;

end