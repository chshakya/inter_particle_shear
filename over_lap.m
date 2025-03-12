function para_overlap=over_lap(d1,d2,r1,r2,gap1,gap2)

% d1=disp_mean(1);
% d2=disp_mean(2);
% r1=particle_dia_upper_perp/2;
% r2=particle_dia_lower_perp/2;
% gap1=gap_mean(1);
% gap2=gap_mean(2);

delta=(gap1-gap2)/1000;
c_2=(((((d2/2)^2)-((d1/2)^2))/delta)+delta)/2;
o_1=(((d1/2)^2+(c_2^2))^0.5)-c_2;
c_calc=c_2+o_1;
c_ass=r1+r2;
error_sum_r=c_calc-c_ass;
o_2=o_1+delta;
para_overlap.error_sum_r=error_sum_r;
para_overlap.o_1=o_1;
para_overlap.o_2=o_2;
para_overlap.sum_radii_calc=c_calc;

end
