% r1=particle_dia_lower_perp/2;
% r2=particle_dia_upper_perp/2;
% disp=Samplegap.rep.processed_data.lin_disp_m(1);
% r=moment_arm;

% function over_lap_r1r2=overlap_max_r1_r2_true(r1,r2,disp)
% over_lap_r1r2=r1+r2-sqrt(((r1+r2)^2)-((disp/2)^2));
% end

 function over_lap_r1r2=overlap_max_r1_r2_true(r1,r2,disp,r)
theta=disp/r;
over_lap_r1r2.d3_noncentric=(r1+r2)-sqrt(((r1+r2)^2)-(((r*sin(theta))/(1+cos(theta)))^2));
d_theo=2*r*tan(disp/(2*r));
over_lap_r1r2.d3_centric=(r1+r2)-sqrt(((r1+r2)^2)-((d_theo/2)^2));
over_lap_r1r2.lin=r1+r2-sqrt(((r1+r2)^2)-((disp/2)^2));
 end