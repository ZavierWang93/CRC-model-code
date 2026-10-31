function funval = cohortObjective2(b1, para)

a0 = para.a0;
a1 = para.a1;
t0 = para.t0;
t1 = para.t1;
b0 = para.b0;
incRate = para.incRate;

numer1 = exp((1/2)*(t0^2*a0+2*t0^2*a1+2*a0*t1*t0+a0*t1^2)*a0/(-t0*a1+t0*a0+a1*t1-a0*t1))*b0*sqrt(pi)*sqrt(2)*erf((1/2)*sqrt(2)*a0/sqrt((-a1+a0)/(t0-t1)))*a1-2*exp((t0*a1+2*a0*t1)*t0*a0/(-t0*a1+t0*a0+a1*t1-a0*t1))*b1*sqrt((-a1+a0)/(t0-t1))-exp((1/2)*(t0^2*a0+2*t0^2*a1+2*a0*t1*t0+a0*t1^2)*a0/(-t0*a1+t0*a0+a1*t1-a0*t1))*b1*sqrt(pi)*sqrt(2)*erf((1/2)*sqrt(2)*a0/sqrt((-a1+a0)/(t0-t1)))*a0+2*exp((t0*a1+2*a0*t1)*t0*a0/(-t0*a1+t0*a0+a1*t1-a0*t1))*b0*sqrt((-a1+a0)/(t0-t1))-exp((1/2)*(t0^2*a0+2*t0^2*a1+2*a0*t1*t0+a0*t1^2)*a0/(-t0*a1+t0*a0+a1*t1-a0*t1))*b0*sqrt(pi)*sqrt(2)*a1*erf((1/2)*sqrt(2)*a1/sqrt((-a1+a0)/(t0-t1)))+2*exp((1/2)*(2*t0^2*a1*a0-t0^2*a1^2+t0^2*a0^2+2*a0^2*t1*t0+2*t0*a1^2*t1+a0^2*t1^2-a1^2*t1^2)/(-t0*a1+t0*a0+a1*t1-a0*t1))*b1*sqrt((-a1+a0)/(t0-t1))+exp((1/2)*(t0^2*a0+2*t0^2*a1+2*a0*t1*t0+a0*t1^2)*a0/(-t0*a1+t0*a0+a1*t1-a0*t1))*b1*sqrt(pi)*sqrt(2)*a0*erf((1/2)*sqrt(2)*a1/sqrt((-a1+a0)/(t0-t1)))-2*exp((1/2)*(2*t0^2*a1*a0-t0^2*a1^2+t0^2*a0^2+2*a0^2*t1*t0+2*t0*a1^2*t1+a0^2*t1^2-a1^2*t1^2)/(-t0*a1+t0*a0+a1*t1-a0*t1))*b0*sqrt((-a1+a0)/(t0-t1));
numer2 = numer1*exp(-(t0*a1+2*a0*t1)*t0*a0/((-a1+a0)*(t0-t1)));
numer = 0.5/(sqrt((a0-a1)/(t0-t1))*(a0-a1))*numer2;

denom = -(1/2)*sqrt(2)*sqrt(pi)*exp((1/2)*(t0-t1)*a0^2./(-a1+a0))*(erf((1/2)*sqrt(2)*a0./sqrt((-a1+a0)/(t0-t1)))-erf((1/2)*sqrt(2)*a1./sqrt((-a1+a0)/(t0-t1))))./sqrt((-a1+a0)/(t0-t1));

funval = (1-(numer./denom/incRate)).^2;
