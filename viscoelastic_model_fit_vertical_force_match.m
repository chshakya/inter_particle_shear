close all 
clear all
cd 'M:\Rheometer\scripts'
% Calculate force on particles using Jian et al, Powder Technology 350
% (2019), p 51.  Equation 14, 17, and 23


%  Based on a standard linear solid (spring in parallel with MAxwell
%  element)
%
% Ge=1.25e+04;  % shear modulus at t->inf
% G1=0.25*Ge; % "relaxable" part of shear modulus 
% eta=1; %viscosity
% vp=0.5;  %Poisson ratio
% T= 100%0.05*(1/(5.6164e-04)); %relaxation time
% R=(18.5/2)/1000;   %sphere radius
% Rstar=R/2;   %1/R=1/R1+1/R2
% str_ain=0.018;
% dmax=str_ain*2*2*R;  % deformation of sphere when they are exactly on top of each other
% y=2*R-dmax;  %vertical distance between centers of spheres


Ge=0.45e+04;  % shear modulus at t->inf
G1=4.0*Ge; % "relaxable" part of shear modulus 
eta=1; %viscosity
vp=0.5;  %Poisson ratio
T= 2000%0.05*(1/(5.6164e-04)); %relaxation time
R=(18.5/2)/1000;   %sphere radius
Rstar=R/2;   %1/R=1/R1+1/R2
str_ain=0.018;
dmax=str_ain*2*2*R;  % deformation of sphere when they are exactly on top of each other
y=2*R-dmax;  %vertical distance between centers of spheres



v_char=sqrt(4*R^2-y^2)/T;
geometric_travel=sqrt(4*R^2-y^2)/T;

v_fact=2.38e-5*[1 2 5 10]%*logspace(-4,4,91);
col=hsv(2*length(v_fact));
mu=5e-3;



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

force_struct(m).x=x;
force_struct(m).Fe=Fe;
force_struct(m).Fe_hz=((x-xo).*Fe./d);

force_struct(m).del=del;

%Dissipation during elastic contact

dissipation_elastic_x(m)=trapz(x,((x-xo).*Fe./d));
 
%Dissipative part:
pref=4*sqrt(Rstar)/(2*(1-vp));
ds=dt/500;  %time step for numerical integration (eq 23)
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

force_struct(m).Fd=Fd;

Res(m,:)=Fe+Fd;
force_struct(m).Res=Fe+Fd;
force_struct(m).Res_x=(x-xo).*(Fe+Fd)./d;
force_struct(m).x=x;

trigger_val=0;
x_zero(m)=x(end);
[max_Res,max_Res_pos]=max(Res(m,:));
[min_Res,min_Res_pos]=min(Res(m,:));
for i=min(min_Res_pos,max_Res_pos)+10:max(min_Res_pos,max_Res_pos)-1
    if or(and(Res(m,i)<0,Res(m,i+1)>0),and(Res(m,i)>0,Res(m,i+1)<0))
        zero_fit=polyfit(x(i:i+1),Res(m,i:i+1),1);
        x_zero_array=x(i:i+1);
        y_zero_array=Res(m,i:i+1);
        x_zero(m)=-zero_fit(2)/zero_fit(1);
        trigger_val=1;
        break
    end
end


if trigger_val==1
    x_cut(m).contact_length=[(x(1:i)),(x_zero(m))];
    res_cut(m).resultant_viscoelastic=[Res(m,1:i),0];
    res_cut(m).resultant_viscoelastic_x=[(x(1:i)-xo).*Res(m,1:i)./d(1:i),0];
else
    x_cut(m).contact_length=x;
    res_cut(m).resultant_viscoelastic=Res(m,:);
    res_cut(m).resultant_viscoelastic_x=(x-xo).*Res(m,:)./d;
end

dissipation(m)=trapz((x_cut(m).contact_length),(res_cut(m).resultant_viscoelastic_x));

%% friction dissipation
res_cut(m).friction=mu*res_cut(m).resultant_viscoelastic;

if trigger_val==1
res_cut(m).friction_x= [mu*y*Res(m,1:i)./d(1:i),0];%y*res_cut(m).friction(m,:)./[d(1:i),0];
else
res_cut(m).friction_x= mu*y*Res(m,:)./d;%y*res_cut(m).friction(m,:)./[d(1:i),0];    
end

res_cut(m)
res_cut(m).resultant_x=res_cut(m).resultant_viscoelastic_x-res_cut(m).friction_x;

dissipation_friction(m)=trapz(x_cut(m).contact_length,res_cut(m).friction_x);


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

%% minimum cycle length
x_geometric(m)=x(end);


%% clear used variables
clear Fe Fd FF x d v dt t del pref ds i s xi di deli d12i ddoti del_t expo integ trigger_val max_Res max_Res_pos ...
    min_Res min_Res_pos zero_fit x_zero_array y_zero_array Array FF5

%% dissipation
% in the direction of contact


end
%% plot vertical forces for paper
f=figure
for i=1:m
    name_tag=sprintf('%0.2e',v_fact(i));
    plot(x_plot(i,:),Resy(i,:)+FyF(i,:),"DisplayName",name_tag,'LineWidth',2)
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
AddLetters2Plots({ax1},{'b)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
saveas(gcf,'M:\paper_writing\inter_particle_shear\ver_force_vs_disp_simulation.svg')
saveas(gcf,'M:\paper_writing\inter_particle_shear\ver_force_vs_disp_simulation','epsc')

%% minimum cycle length
normalized_x=x_zero./x_geometric;
normalized_vel=v_fact/v_char;
[x_min,x_min_pos]=min(normalized_x);
vel_min_x_norm=normalized_vel(x_min_pos);
% 
% 
%%
figure
scatter(v_fact,x_zero)

%%
figure
plot(Resx(1,:))
hold on
plot(Resx(1,:)-FxF(1,:))
yline(0)
hold off

figure
plot(Resy(1,:))
hold on
plot(Resy(1,:)+FyF(1,:))
yline(0)
hold off


%% plots

% figure
% plot(x_plot(m,:),force_struct(m).Fe,'DisplayName','Hertz Model')
% hold on
% for m =1:length(v_fact)
% %   plot(x_plot(m,:),Res(m,:),'DisplayName', strcat('velocity ',string(compose("%5.2e",v*v_fact(m)))),'LineStyle','--')
%     plot(x_plot(m,:),Res(m,:),'DisplayName', strcat('velocity ',string(compose("%5.2e",v_fact(m)))),'LineStyle','-','Color',col(m+1,:))  
%     
%     hold on
%     plot(x_zero_array,y_zero_array,'o')
%     plot(x_plot(m,:),Res_mod(m,:),'DisplayName', strcat('velocity ',string(compose("%5.2e",v_fact(m)))),'LineStyle','--','Color',col(m,:))
% 
%   
%   
% %   xline(x_end(m),'Color',col(m,:))
% end
% % yline(0)
% hold off
% 
% figure 
% plot(x_plot(m,:),FxH(m,:),'DisplayName','Hertz Model')
% hold on
% for m =1:length(v_fact)
% % plot(x_plot(m,:),FxD(m,:),'DisplayName',strcat('distance scale  ',string(compose("%5.2e",v*v_fact(m)*T))),'Color',col(m,:))
% plot(x_plot(m,:),Resx(m,:)-FxF(m,:),'DisplayName', strcat('velocity ',string(compose("%5.2e",v*v_fact(m)))),'LineStyle','--')
% 
% 
% end
% yline(0)
% xline(0)
% hold off
% xlabel("Displacement")
% ylabel("Horizontal component of forces")
% legend('Location','southeast')
% set(gca,'FontSize',12)
% 
% saveas(gcf,'/home/ubuntu_pcc3/Downloads/matlab_script-examples/plots/dissipative','epsc')
 
% 
% 
% figure 
% plot(x_plot(m,:),FyH(m,:),'DisplayName','Hertz Model')
% hold on
% for m =1:length(v_fact)
% % plot(x_plot(m,:),FxD(m,:),'DisplayName',strcat('distance scale  ',string(compose("%5.2e",v*v_fact(m)*T))),'Color',col(m,:))
% plot(x_plot(m,:),Resy(m,:)-FyF(m,:),'DisplayName', strcat('velocity ',string(compose("%5.2e",v*v_fact(m)))),'LineStyle','--')
% 
% 
% end
% yline(0)
% xline(0)
% hold off
% xlabel("Displacement")
% ylabel("Vertical component of forces")
% legend('Location','southeast')
% %figure,hold on
% %plot(x,FxH,'b')
% %plot(x,FxD+FxH,'r')

% x vs. Hz elastic force
figure
plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Elastic')
yline(0,'HandleVisibility','off')
xline(force_struct(1).x(end)/2,'HandleVisibility','off')
hold off
legend
xlabel('displacement (m)')
ylabel('Forces in Horizontal direction(N)')
% ylabel({'Elastic'; 'normalized';'by';'maximum(1)'},'HorizontalAlignment','right')
% set(get(gca,'ylabel'),'rotation',0)

% x vs. Hz force
figure
plot(force_struct(1).x,force_struct(1).del/max(force_struct(1).del),'DisplayName','Geometeric overlap')
hold on
plot(force_struct(1).x,force_struct(1).Res/max(force_struct(1).Res),'DisplayName','Forces Viscous+elastic')
yline(0,'HandleVisibility','off')
hold off
legend
xlabel('dispacement (m)')
ylabel({'Value'; 'normalized';'by';'maximum(1)'},'HorizontalAlignment','right')
set(get(gca,'ylabel'),'rotation',0)




% %% contact length variation discrete
% figure
% loglog(v_fact,x_end,'DisplayName','Simulation Data')
% hold on
% xline(v_char,'DisplayName','Chracteristic velocity')
% hold off
% xlabel('velocity of moving probe(m/s)')
% ylabel('Length of contact(m)')
% title('Variation of contact length for 18.5 \phi spheres inter-particle shear')
% legend

%% normalized
figure
loglog(normalized_vel,normalized_x,'DisplayName','Simulation Data')
hold off
xlabel('velocity of moving probe/characteristic velocity')
ylabel({'Resultant Length'; 'of contact/'; 'geometric';'contact length'},'HorizontalAlignment','right')
set(get(gca,'ylabel'),'rotation',0)
title('Variation of contact length for 18.5 \phi spheres inter-particle shear')


figure
loglog(normalized_vel,dissipation/max(abs(dissipation)),'DisplayName','Dissipation')
yyaxis left
xlabel('velocity of moving probe/characteristic velocity')
ylabel('Dissipation x(J)/max dissipation(J)')
hold on
loglog(normalized_vel,normalized_x,'DisplayName','normailzed contact length')
yyaxis right
hold on
xline(vel_min_x_norm,'DisplayName','minimum contact length')
hold off
ylabel('Resultant Length of contact/geometric contact length')
set(get(gca,'ylabel'),'rotation',0)
title('contact length variation for 18.5mm \phi spheres')
legend

figure
loglog(normalized_vel,abs(dissipation),'DisplayName','Viscous')
hold on
loglog(normalized_vel,abs(dissipation_friction),'DisplayName','Friction')
loglog(normalized_vel,abs(dissipation_friction)+abs(dissipation),'DisplayName','Effective')
hold off
xlabel('velocity of moving probe/characteristic velocity')
ylabel('Dissipation x(J)')
set(get(gca,'ylabel'),'rotation',0)

figure
loglog(normalized_vel,abs(dissipation_friction),'DisplayName','Friction')
xlabel('velocity of moving probe/characteristic velocity')
ylabel('Dissipation x(J)')
set(get(gca,'ylabel'),'rotation',0)

figure
loglog(normalized_vel,abs(dissipation)./abs(dissipation_friction),'DisplayName','Viscous/Friction')
hold on
yline(1)
xline(1)
xline(10)
hold off
xlabel('velocity of moving probe/characteristic velocity')
ylabel('Viscous Dissipation x/Friction dissipation x')
set(get(gca,'ylabel'),'rotation',0)

figure
plot(normalized_x(1:floor(end/2)),abs(dissipation(1:floor(end/2))))
hold on
plot(normalized_x(floor(end/2): end),abs(dissipation(floor(end/2): end)))
hold off
set(get(gca,'ylabel'),'rotation',0)

figure
loglog(normalized_vel,abs(dissipation))


%%% plots for presentation

% figure;plot(x,-FxH(1,:));hold on;yline(0);xline(x(end)/2);hold off;xlabel('Displacement (m)');ylabel('Hz Forces (N)')

figure
plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Hertz');hold on;
plot(force_struct(1).x,-force_struct(1).Res_x,'DisplayName','Dissipative+Hertz')
yline(0,'HandleVisibility','off')
hold off
xlabel('Displacement (m)')
ylabel('Horizontal Force (N)')

figure
plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Hertz');hold on;
plot(x_cut(1).contact_length,-res_cut(1).resultant_viscoelastic_x,'DisplayName','Dissipative+Hertz')
yline(0,'HandleVisibility','off')
hold off
xlabel('Displacement (m)')
ylabel('Horizontal Force (N)')

figure
plot(force_struct(1).x/R,-force_struct(1).Fe_hz/dmax,'DisplayName','Hertz');hold on;
plot(x_cut(1).contact_length/R,-res_cut(1).resultant_x/dmax,'DisplayName','Hertz+Friction')
plot(x_cut(1).contact_length/R,-res_cut(1).friction_x/dmax,'DisplayName','Friction')
yline(0,'HandleVisibility','off')
xline(force_struct(1).x(end)/2,'HandleVisibility','off')
hold off
xlabel('^{Displacement}/_{Avg. Radius}')
ylabel('^{Horizontal Force(N)}/_{Max. overlap(m)}')
legend


figure
plot(force_struct(1).x/R,-force_struct(1).Fe_hz/dmax,'DisplayName','Hertz','LineStyle','--');
hold on;
for i=1:length(v_fact)
plot(x_cut(i).contact_length/R,-res_cut(i).resultant_x/dmax,'HandleVisibility','off')
end
yline(0,'HandleVisibility','off')
hold off
xlabel('^{Displacement}/_{Avg. Radius}')
ylabel('^{Horizontal Force(N)}/_{Max. overlap(m)}')

figure
plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Hertz');hold on;
plot(x_cut(1).contact_length,-res_cut(1).resultant_x,'DisplayName','Hertz+Dissipation+Friction')
plot(x_cut(1).contact_length,-res_cut(1).friction_x,'DisplayName','Friction')
yline(0,'HandleVisibility','off')
xline(force_struct(1).x(end)/2,'HandleVisibility','off')
hold off
xlabel('Displacement (m)')
ylabel('Horizontal Force(N)')
legend('Location','southwest')


figure
plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Hertz');hold on;

plot(x_cut(1).contact_length,-res_cut(1).resultant_viscoelastic_x,'DisplayName','Dissipation+Hertz')
yline(0,'HandleVisibility','off')
xline(force_struct(1).x(end)/2,'HandleVisibility','off')
hold off
xlabel('Displacement (m)')
ylabel('Horizontal Force (N)')
legend('Location','southwest')
%% plots for paper
figure
subplot(2,2,1)
plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Hertz')
[ymin,idx_min] = min(-force_struct(1).Fe_hz) ;
[ymax,idx_max] = max(-force_struct(1).Fe_hz) ; 
text(force_struct(1).x(idx_min),ymin,['t: ' num2str(ymin)]);
text(force_struct(1).x(idx_max),ymax,['p: ' num2str(ymax)]);

yline(0,'HandleVisibility','off')
xline(force_struct(1).x(end)/2,'HandleVisibility','off')
hold off
legend('Location','southwest')
xlabel('displacement (m)')
ylabel('Horizontal Force(N)')

subplot(2,2,2)
plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Hertz');hold on;
plot(x_cut(1).contact_length,-res_cut(1).resultant_viscoelastic_x,'DisplayName','Dissipation+Hertz')
yline(0,'HandleVisibility','off')
xline(force_struct(1).x(end)/2,'HandleVisibility','off')
hold off
xlabel('Displacement (m)')
ylabel('Horizontal Force (N)')
legend('Location','southwest')

subplot(2,2,3)
plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Hertz');hold on;
plot(x_cut(1).contact_length,-res_cut(1).resultant_x,'DisplayName','Hertz+Dissipation+Friction')
plot(x_cut(1).contact_length,-res_cut(1).friction_x,'DisplayName','Friction')
yline(0,'HandleVisibility','off')
xline(force_struct(1).x(end)/2,'HandleVisibility','off')
hold off
xlabel('Displacement (m)')
ylabel('Horizontal Force(N)')
legend('Location','southwest')


subplot(2,2,4)
% plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Hertz','LineStyle','--');
% hold on;
for i=1:length(v_fact)
    name_tag=sprintf('%d',v_fact(i))
plot(x_cut(i).contact_length,-res_cut(i).resultant_x,'DisplayName',name_tag)
hold on
end
yline(0,'HandleVisibility','off')
hold off
xlabel('Displacement (m)')
ylabel('Horizontal Force(N)')
lgd = legend('Location','southwest');
lgd.Title.String = 'Rotational Speed (m/s)';

%% Hz force simulation
f=figure
plot(force_struct(1).x,-force_struct(1).Fe_hz,'DisplayName','Hertz');
hold on;
plot(x_cut(1).contact_length,-res_cut(1).resultant_x,'DisplayName','Hertz+Dissipation+Friction')
yline(0,'HandleVisibility','off')
xline(force_struct(1).x(end)/2,'HandleVisibility','off')
hold off
[ymin,idx_min] = min(-res_cut(1).resultant_x) ;
[ymax,idx_max] = max(-res_cut(1).resultant_x) ; 
text(force_struct(1).x(idx_min),ymin,['t: ' num2str(round(ymin,4))],'FontSize',18);
text(force_struct(1).x(idx_max),ymax,['p: ' num2str(round(ymax,4))],'FontSize',18);
% xlim([0.0045 0.0175])
xlabel("$\theta$ [m]",'Interpreter','latex')
ylabel("$F_h$ [N]",'Interpreter','latex')
ax1 = gca;
ax1.LineWidth=2;
ax1.FontSize = 20;
ax1.YAxis.Exponent = -3;
ax1.XAxis.Exponent = -3;
axis('square')
f.Position=[200 300 560 500]
AddLetters2Plots({ax1},{'a)'}, 'HShift', 0.06, 'VShift', 0, 'FontSize', 18)
% saveas(gcf,'M:\paper_writing\inter_particle_shear\hz_force_vs_disp_ampdiff.svg')
% saveas(gcf,'M:\paper_writing\inter_particle_shear\hz_force_vs_disp_ampdiff','epsc')

