function cdf_v = cdfUnderlyingLogFIT(vec_v, stage, contractor_idx, model_specs, age)

aa_base=model_specs.mu{contractor_idx}(stage);
age_eff=model_specs.age_eff_mean{contractor_idx}(stage);
age_delta = max(age - 50, 0);
aa = aa_base + age_eff * age_delta;
aa = max(aa, 1e-6);

bb=model_specs.sigma{contractor_idx}(stage);
[cc, dd] = lognormalParameters(aa, bb);
pZero=model_specs.pZero{contractor_idx}(stage);

cdf_v = cdf('Lognormal', vec_v, cc, dd);
cdf_v = pZero+(1-pZero)*cdf_v;
