function cdf_h = cdfObservedLogFIT(vec_h, stage, contractor_idx, model_specs, age)

omega_eps = 1e-16;

numPts = 20;
next_iter = true;
convergence_tol = 1e-3;
old_cdf_h = 1;
while (next_iter)
    [vec_v, pdf_v, dv] = pdfNonzeroLogFIT(stage,model_specs,numPts,contractor_idx,age);

    vec_omega = getOmega(vec_v, model_specs, contractor_idx);
    q0 = model_specs.q0{contractor_idx};
    pZero = model_specs.pZero{contractor_idx}(stage);

    cdf_h_given_v = zeros(numel(vec_h), numel(vec_v));
    for ii=1:numel(vec_v)
        val_v = vec_v(ii);
        val_q0 = q0;
        val_omega = vec_omega(ii);
        if (val_omega<omega_eps)
            zz = zeros(numel(vec_h),1);
        else

            zz = 1-0.5*(erfc((vec_h-val_v)./(val_omega*sqrt(2))));
        end
        cdf_h_given_v(:,ii) = pZero+(1-pZero)*(val_q0+(1-val_q0)*zz);
    end
    cdf_h = cdf_h_given_v*pdf_v*dv;

    int_diff = abs(1-cdf_h./old_cdf_h);
    if (all(int_diff<convergence_tol))
        next_iter = false;
    else
        numPts = 2*numPts;
        old_cdf_h = cdf_h;
    end
end
