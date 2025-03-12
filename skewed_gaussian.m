function gau=skewed_gaussian(x,y,sd)

[max_val,max_pos]=max(y);
init(2)=x(max_pos); %center posn of peak
init(1)=sd; %stand deviation , scale
init(3)=0.1; %right skewed a>0, left skewed a<0
init(4)=max_val;

tbl = table(x', y');
modelfun =@(A,x) A(4)*(2/(A(1)*sqrt(pi()))).*exp(-((x-A(2)).^2)/(2*(A(1)^2))).*(0.5*(1+erf((A(3)*(x-A(2))/A(1))/sqrt(2))));
beta0 = [init(1); init(2); init(3); init(4)]; %parameter fit range
mdl = fitnlm(tbl, modelfun, beta0);
A=mdl.Coefficients.Estimate;
func=A(4)*(2/(A(1)*sqrt(pi()))).*exp(-((x-A(2)).^2)/(2*(A(1)^2))).*(0.5*(1+erf((A(3)*(x-A(2))/A(1))/sqrt(2))));
gau=mdl;
figure
scatter(x,y)
hold on
plot(x,func)
end



