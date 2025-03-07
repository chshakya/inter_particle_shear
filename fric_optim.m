function optim_fric=fric_optim(mu1,mu2,x_zero_exp,x,hz_comp,norm_force,cosine,divisions)
%% friction coeffcient optimization

div=1:divisions;
mu_val=mu1+((mu2-mu1)*div/max(div));
for j=1:max(div)
    Fric_tot=mu_val(j)*norm_force;
    Fric_hz=Fric_tot.*cosine;
    Eff_hz=hz_comp-Fric_hz;
    [max_Hz,max_Hz_pos]=max(Eff_hz);
    [min_Hz,min_Hz_pos]=min(Eff_hz);
    min(min_Hz_pos,max_Hz_pos)
    max(min_Hz_pos,max_Hz_pos)
        for i=min(min_Hz_pos,max_Hz_pos):max(min_Hz_pos,max_Hz_pos)
            if and(Eff_hz(i)<0,Eff_hz(i+1)>0)
                flip_index=i;
                zero_fit=polyfit(x(i:i+1),Eff_hz(i:i+1),1);
                x_zero_num=-zero_fit(2)/zero_fit(1);
                break
            end
        end
     diff_x_zero(j)= x_zero_exp-x_zero_num; 
end

% figure
% plot(diff_x_zero,mu_val)
% xlabel("Distance of experimental flip and numeric flip (m)")
% ylabel("\mu")
mu_optim_fit=polyfit(diff_x_zero,mu_val,1)
optim_fric=mu_optim_fit(2);

end
