function [mu, sigma] = lognormalParameters(m, s)

mu = log(m.^2./sqrt(s.^2+m.^2));
sigma = sqrt(log(s.^2./m.^2+1));

[mm, vv] = lognstat(mu, sigma);

eps2 = 1e-6;
VV = (mm-m).^2+(sqrt(vv)-s).^2;
Vmax = max(VV);
if (Vmax>eps2)
    error(['Error in lognormalParameters: Vmax = ' num2str(Vmax)]);
end;