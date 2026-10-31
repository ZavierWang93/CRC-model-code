function  [vec_v, pdf_v, dv]=pdfNonzeroLogFIT(stages,paras,numpts,lab,age)

  aa_base=paras.mu{lab}(stages);
  age_eff=paras.age_eff_mean{lab}(stages);
    age_delta = max(age - 50, 0);
    aa = aa_base + age_eff * age_delta;
  aa = max(aa, 1e-6);

  bb=paras.sigma{lab}(stages);
  [cc, dd] = lognormalParameters(aa, bb);

  eps = 1e-100;
  cdf_err = 1e-10;

  min_v = icdf('Lognormal', cdf_err, cc, dd);
  max_v = icdf('Lognormal', 1-cdf_err, cc, dd);

    min_v = max(eps,min_v);
    max_v = min(log(1e12),max_v);
    max_v = max(max_v, min_v+eps);

    vec_v = linspace(min_v, max_v, numpts);

    pdf_v = pdf('Lognormal', vec_v, cc, dd)';
    pdf_v = pdf_v+eps;

    dv = vec_v(2)-vec_v(1);

    normC = sum(pdf_v)*dv;

    if (normC>0.99)
        pdf_v = pdf_v/normC;
    end

    vec_v = vec_v';

    if (any(isnan(vec_v)) || any(isnan(pdf_v)))
        error('Error in pdf_logFIT_underlying_nonzero');
    end

end
