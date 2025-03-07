function over_lap_r1r2=overlap_r1_r2_true(r1,r2,disp)
over_lap_r1r2=r1+r2-sqrt(((r1+r2)^2)-((disp/2)^2));
end