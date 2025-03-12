function para=skewed_parabola(x,y,sd)

[max_val,max_pos]=max(y);

fit_poly=polyfit(x,y,2);
init(1)=fit_poly(1); %amplitude A(1)
init(2)=(x(1)-x(end))/2; %center posn of peak x dirction A(2)
init(3)=0.1; %slope of line
init(4)=max_val;
init(5)=y(1);

tbl = table(x', y');
modelfun =@(A,x) A(1)*((x-A(2)).^2) + x*A(3) + A(4)+ A(5);
beta0 = [init(1); init(2); init(3); init(4); init(5)]; %parameter fit range
mdl = fitnlm(tbl, modelfun, beta0);
A=mdl.Coefficients.Estimate;
func= A(1)*((x-A(2)).^2) + x*A(3) + A(4)+ A(5);
para=mdl;
figure
scatter(x,y)
hold on
plot(x,func)
end



