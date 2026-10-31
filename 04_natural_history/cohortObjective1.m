function funval = cohortObjective1(a1, para)

a0 = para.a0;
t0 = para.t0;
t1 = para.t1;
mRate = para.mRate;

numer = 1-exp(-(1/2)*(t1-t0)*(-a1*t1-a0*t1+t0*a1+t0*a0)/(t0-t1));
denom = -(1/2)*sqrt(2)*sqrt(pi)*exp((1/2)*(t0-t1)*a0^2./(-a1+a0)).*(erf((1/2)*sqrt(2)*a0./sqrt((-a1+a0)/(t0-t1)))-erf((1/2)*sqrt(2)*a1./sqrt((-a1+a0)/(t0-t1))))./sqrt((-a1+a0)/(t0-t1));

funval = (mRate-numer./denom).^2;
