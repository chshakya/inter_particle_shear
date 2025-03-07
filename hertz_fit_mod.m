function hertz_friction=hertz_fit_mod(x,y,ver,time_array,v_overlap,r1,r2,overlap_exp,E,nu,mu1,mu2,x_zero_exp,v_hz_mean,max_lin_disp_m)
% x: displacement [m].
% y: Horizontal force [N].
% ver: Vertical force [N].

mu=mu1;
%equivalent r
r=(r1*r2)/(r1+r2);

% disp_x=abs(x(1)-x(end));
disp_x=max_lin_disp_m;

% a_exp=disp_x*r2/(2*(r1+r2));
overlap_num=(r1+r2)-(sqrt(((r1+r2)^2)-((disp_x/2)^2)));
hertz_friction.correction_overlap=overlap_exp-overlap_num;
overlap_exp=overlap_num;
d_c_c=r1+r2-overlap_exp; %distance between particle centers

x_m=squeeze(x); %position of moving particles
y_m=d_c_c;

x_f=disp_x/2;
% (x_m(1)+x_m(end))/2; % position of stationary particles
y_f=0;

%for calculation of trig ratios
perp=(x_m-x_f);
base=(y_m-y_f);

%dist between particle centers
dist=sqrt((x_m-x_f).^2+(y_m-y_f).^2);

%trig ratios
sine=perp./dist;
cosine=base./dist;
tangent=perp./base;
alpha=atan(tangent);


%maximum displacement
max_dist=2*sqrt((2*(r1+r2)*overlap_exp)+overlap_exp^2);

%calculating possible center to center distances
dist_pos=dist.*(dist<(r1+r2));
overlap=(r1+r2)-dist_pos; %calculating overlap
a=sqrt(r*overlap); % calculating length of contact if gels were perpendicular

%%elastic constants

%% conversion to shear moduli
G=E/(2*(1+nu));
Ge=G*0.78913;


E_asterisk=E/(2*(1-(nu^2)));


%force calculations
F_elastic_total=(4*Ge*sqrt(r)/(3*(1-nu)))*overlap.^(3/2).*(overlap<max(r1,r2));

F_elastic_vert=F_elastic_total.*cosine;
F_elastic_hor=F_elastic_total.*sine;

% figure
% plot(x_m,dist, "DisplayName" ,"Distance Between particle centers")
% hold on
% plot(x_m,dist_pos, "DisplayName" ,"Distance Between particle centers that are possible")
% plot(x_m,overlap,"DisplayName" ,"Overlaps possible")
% hold off
% xlabel("Distance travelled in by moving particle in x axis")
% ylabel("Relevant Distance")

peak_mismatch=max(F_elastic_hor)+min(F_elastic_hor);

% figure
% plot(x_m,F_elastic_vert,"DisplayName","Vertical component of elastic forces")
% hold on
% plot(x_m, F_elastic_hor,"DisplayName","Horizontal component of elastic forces")
% plot(x_m, F_elastic_total,"DisplayName","Resultant elastic forces")
% yline(0)
% xline(x_f)
% hold off
% xlabel("Distance travelled in x axis")
% ylabel("Forces predicted by herzian model")
% legend

%%considering friction
F_friction_total=mu*F_elastic_total;
F_friction_vert=F_friction_total.*sine;
F_friction_hor=F_friction_total.*cosine;

%effective
F_hor_eff=F_elastic_hor-F_friction_hor;
F_vert_eff=F_elastic_vert-F_friction_vert;

% figure
% plot(x_m, F_elastic_hor,"DisplayName","Hertzian")
% hold on
% yline(0)
% xline(x_f)
% hold off
% xlabel("Distance travelled in x axis")
% ylabel("Horizontal Component of Forces")
% legend('Location','northwest')
% set(gca,'FontSize',12)
% saveas(gcf,'/home/ubuntu_pcc3/Downloads/matlab_script-examples/plots/hertz','epsc')

% 
% figure
% plot(x,y,"DisplayName","Experimental")
% hold on
% plot(x_m, F_elastic_hor,"DisplayName","Hertzian")
% 
% plot(x_m, F_friction_hor,"DisplayName","Frictional")
% plot(x_m, F_hor_eff,"DisplayName","Resultant")
% yline(0)
% xline(x_f)
% hold off
% xlabel("Distance travelled in x axis")
% ylabel("Horizontal Component of Forces")
% legend('Location','northwest')
% 

% set(gca,'FontSize',12)
% saveas(gcf,'/home/ubuntu_pcc3/Downloads/matlab_script-examples/plots/frictional','epsc')


% figure
% plot(x,ver,"DisplayName","Experimental")
% hold on
% plot(x_m, F_elastic_vert,"DisplayName","Herzian")
% plot(x_m, F_friction_vert,"DisplayName","Friction")
% plot(x_m, F_vert_eff,"DisplayName","Resultant")
% yline(0)
% xline(x_f)
% hold off
% xlabel("Distance travelled in x axis (m)")
% ylabel("Vertical component of forces (N)")
% legend

%finding location of displacements
x_contact_index=find(dist_pos,1,'first');
x_detach_index=find(dist_pos,1,'last');
x_contact=x_m(x_contact_index);
x_detach=x_m(x_detach_index);
x_first_quartile_index=x_contact_index+ceil((x_detach_index-x_contact_index)/4);
x_third_quartile_index=x_contact_index+ceil((x_detach_index-x_contact_index)*(3/4));
[max_hz_force_elastic,max_hz_force_elastic_index]=max(F_elastic_hor);
[min_hz_force_elastic,min_hz_force_elastic_index]=min(F_elastic_hor);

%% finding cross point

[max_dF_dx,max_dF_dx_pos]=max(F_hor_eff);
[min_dF_dx,min_dF_dx_pos]=min(F_hor_eff);
for i=min(min_dF_dx_pos,max_dF_dx_pos):max(min_dF_dx_pos,max_dF_dx_pos)
    if and(F_hor_eff(i)<0,F_hor_eff(i+1)>0)
        flip_index=i;
        zero_fit=polyfit(x_m(i:i+1),F_hor_eff(i:i+1),1);
        x_zero=-zero_fit(2)/zero_fit(1);
        break
    end
end
% figure
% plot(x_m, F_hor_eff,"DisplayName","Resultant")
% hold on
% 
% plot(x_zero,0,'r*')
% yline(0)
% xline(x_f)
% hold off
% legend


%% zero cross misfit
hertz_friction.hertz.zero_mismatch=x_f-x_zero_exp; %mismatch in crossing point

%% residual of numeric and experimental data
hor_resi=y-F_hor_eff;
% figure
% plot(x,hor_resi)

ver_resi=ver-F_vert_eff;
% figure
% plot(x,ver_resi)


%ratio theoretical to ecp
% figure
% plot(x(30:end-30),y(30:end-30)./F_hor_eff(30:end-30))
 error_hz=mean(y(30:end-30)./F_hor_eff(30:end-30));
% hold on
% plot(x(30:end-30),ver(30:end-30)./F_vert_eff(30:end-30))
 error_vert=mean(ver(30:end-30)./F_vert_eff(30:end-30));
% 
% hold off
% title("Error Trend")


%% viscoelastic


G1=(G-Ge);
T= 0.01*(1/(5.6164e-04)); %relaxation time
v=


dmax=overlap_exp;  % deformation of sphere when they are exactly on top of each other
y=d_c_c;  %vertical distance between centers of spheres
mu=0.5e-3;

to=time_array(end)-time(1);%(max_dist)/v_hz_mean; % time for crossing
no_datpoints=length(x);

t=time_array; 
x=x_m; %x-position of sphere
xo=x_f;%x-position of fixed sphere
d=dist_pos; 
del=overlap; %indentation

dt=abs(avg(diff(time_array)));
 

%% Dissipative part:
pref=4*sqrt(r)/(3*(1-nu));
ds=dt/50;  %time step for numerical integration (eq 23)
Fd=zeros(size(t));
 
for i=2:length(t)
    s=[0:ds:t(i)];
    xi=v_overlap(i)*(s);
    di=sqrt(y^2+(xi-xo).^2);
    deli=(r1+r2-di);
    d12i=deli.^(1/2);
    ddoti=-(xi-xo)*v_overlap(i)./(2*di);
    del_t=t(i)-s;
    expo=exp(-(del_t)/T);
    integ(i).arr=G1*exp(-(del_t)/T).*d12i.*ddoti;
    Fd(i)=pref*trapz(s,integ(i).arr); %numerical integral for Fd
end

Res(m,:)=F_elastic_total+Fd;

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






%% optimization mu
divisions=100;
mu_optim=fric_optim(mu1,mu2,x_zero_exp,x,F_elastic_hor,F_elastic_total,cosine,divisions);

%%considering friction
F_friction_total=mu_optim*F_elastic_total;
F_friction_vert=F_friction_total.*sine;
F_friction_hor=F_friction_total.*cosine;

%effective
F_hor_eff=F_elastic_hor-F_friction_hor;
F_vert_eff=F_elastic_vert-F_friction_vert;

figure
plot(x,y,"DisplayName","Experimental")
hold on
plot(x_m, F_elastic_hor,"DisplayName","Hertzian")

plot(x_m, F_friction_hor,"DisplayName","Frictional")
plot(x_m, F_hor_eff,"DisplayName","Resultant")
yline(0)
xline(x_f)
hold off
xlabel("Distance travelled in x axis")
ylabel("Horizontal Component of Forces")
title("Optimized Friction coefficient")
legend('Location','northwest')


%% residual of optimized numeric and experimental data
optim_hor_resi=y-F_hor_eff;
% figure
% plot(x,hor_resi)

optim_ver_resi=ver-F_vert_eff;
% figure
% plot(x,ver_resi)


%% mu from experimental data
mu_exp=((y.*cosine)+(ver.*sine))./((y.*sine)+(ver.*cosine));
% figure
% plot(x_m,mu_exp)
% title('Mu distribution')
mu_hard=y./ver;


%% error from optimized mu
error_hz_optim=mean(y(30:end-30)./F_hor_eff(30:end-30));
error_vert_optim=mean(ver(30:end-30)./F_vert_eff(30:end-30));


%% return variables
hertz_friction.hertz.hz=F_elastic_hor;
hertz_friction.hertz.ver=F_elastic_vert;
hertz_friction.hertz.tot=F_elastic_total;
hertz_friction.hertz.sine=sine;
hertz_friction.hertz.cosine=cosine;
hertz_friction.fric.hz=F_friction_hor;
hertz_friction.fric.ver=F_friction_vert;
hertz_friction.fric.eff_hz=F_hor_eff;
hertz_friction.fric.eff_ver=F_hor_eff;
hertz_friction.fric.hz_resi=hor_resi;
hertz_friction.fric.ver_resi=ver_resi;
hertz_friction.overlap_num=overlap_num;
hertz_friction.mu_optim=mu_optim;
hertz_friction.error_hz=error_hz;
hertz_friction.error_ver=error_vert;
hertz_friction.optim.hz_resi=optim_hor_resi;
hertz_friction.optim.ver_resi=optim_ver_resi;
hertz_friction.optim.error_hz=error_hz_optim;
hertz_friction.optim.error_ver=error_vert_optim;
hertz_friction.exp.mu_exp=mu_exp;
hertz_friction.exp.mu_hard=mu_hard;
hertz_friction.geom.alpha=alpha;
hertz_friction.geom.tangent=tangent;
end
