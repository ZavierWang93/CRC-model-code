function  funval=incidenceConstraint(model_specs, mc)

ageVec = mc.knot_incidence_rate(1):mc.knot_incidence_rate(end);

incidence_rate(1, :) = model_specs.inc_NA;
incidence_rate(2, :) = model_specs.inc_LSL;
AA{1}=  pchip(mc.knot_incidence_rate, incidence_rate(1, :), ageVec);
AA{2}=  pchip(mc.knot_incidence_rate, incidence_rate(2, :), ageVec);

df1 = diff(AA{1});
df2 = diff(AA{2});

is_monotonic1 = all(df1 > 0);
is_monotonic2 = all(df2 > 0);

if is_monotonic1~=1 || is_monotonic2~=1
    funval=-100000;
else
    funval=0;
end

end