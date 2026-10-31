function output = integrateFITLikelihood(data_this,stages,model_specs,mc,lab,varargin)

    ZERO=mc.ZERO;
    VAILD=mc.EXACT;
    INVAILD=mc.INVALID;
    R_CENSORED=mc.RIGHT_CENSORED;
    L_CENSORED=mc.LEFT_CENSORED;

    target_vec = [ZERO VAILD R_CENSORED L_CENSORED INVAILD];
    target = target_vec;

    log_L_censored=mc.logFIT_left_censored(lab);
    log_R_censored=mc.logFIT_right_censored(lab);

    knot_omega=model_specs.knot_omega;
    omega=model_specs.omega;
    pZero=model_specs.pZero(stages);
    qNA=model_specs.qNA;
    q0=model_specs.q0;

    C_pZero = 1 - pZero; C_q0 = 1 - q0; C_qNA = 1 - qNA;
    knot_x = log(1 + knot_omega);
    pp_omega = pchip(knot_x, omega);

    mu_50 = model_specs.mu(stages);
    age_eff = model_specs.age_eff_mean(stages);
    bb = model_specs.sigma(stages);

    log_NA = log(1 + mc.INVALID_FIT);
    inv_SQRT2 = 0.7071067811865476;
    inv_SQRT2PI = 0.3989422804014327;

    rel_tol = 1e-6;
    if isfield(mc, 'riemann_convergence_tol')
        rel_tol = mc.riemann_convergence_tol;
    end

    prob_individual = [];
    age_individual = [];
    id_individual = [];
    row_order_individual = [];
    CRC_stage_individual = [];
    last_index = 0;

    for target_1 = target
        for target_2 = target
            data_this_target = data_this{target_1, target_2};
            if isempty(data_this_target), continue; end

            h1_vec = data_this_target(:, 1);
            h2_vec = data_this_target(:, 2);
            age_vec = data_this_target(:, 3);

            is_na1 = (h1_vec(1) == log_NA);
            is_na2 = (h2_vec(1) == log_NA);

            n_obs = size(data_this_target, 1);

            age_delta = max(age_vec - 50, 0);
            aa_vec = mu_50 + age_eff * age_delta;
            aa_vec = max(aa_vec, 1e-6);

            cc_vec = log(aa_vec.^2 ./ sqrt(bb.^2 + aa_vec.^2));
            dd_vec = sqrt(log(bb.^2 ./ aa_vec.^2 + 1));

            two_dd_sq_inv_vec = 1 ./ (2 * dd_vec.^2);
            dd_sqrt2pi_inv_vec = 1 ./ (dd_vec * 2.506628274631000);

            cdf_err = 1e-12;
            min_v_all = zeros(n_obs, 1);
            max_v_all = zeros(n_obs, 1);
            for ii = 1:n_obs
                min_v_all(ii) = max(1e-100, icdf('Lognormal', cdf_err, cc_vec(ii), dd_vec(ii)));
                max_v_all(ii) = min(log(1e7), icdf('Lognormal', 1-cdf_err, cc_vec(ii), dd_vec(ii)));
            end
            min_v = min(min_v_all);
            max_v = max(max_v_all);

            integrand = @(v) computeIntegrandCore(v, h1_vec, h2_vec, target_1, target_2, ...
                cc_vec, dd_sqrt2pi_inv_vec, two_dd_sq_inv_vec, ...
                knot_x, omega, pp_omega, ...
                pZero, q0, qNA, C_pZero, C_q0, C_qNA, ...
                VAILD, ZERO, R_CENSORED, L_CENSORED, INVAILD, ...
                is_na1, is_na2, inv_SQRT2, inv_SQRT2PI);

            prob_this = integral(integrand, min_v, max_v, ...
                'ArrayValued', true, ...
                'AbsTol', 1e-10, 'RelTol', rel_tol);

            vv = last_index + (1:numel(prob_this));
            prob_individual(vv) = max(prob_this, 1e-100);
            age_individual(vv) = data_this_target(:,3);
            id_individual(vv) = data_this_target(:,5);

            if (stages == mc.CRC)
                CRC_stage_individual(vv) = data_this_target(:,6);
            end

            if size(data_this_target, 2) >= 7
                row_order_individual(vv) = data_this_target(:,7);
            end

            last_index = last_index + numel(prob_this);
        end
    end

    output.prob = prob_individual';
    output.age = age_individual';
    output.screenee_id = id_individual';

    if (stages == mc.CRC)
        output.CRC_stage = CRC_stage_individual';
    end

    if ~isempty(row_order_individual)
        output.row_order = row_order_individual';
    end

    output.final_numpts_mat = zeros(5,5);
end

function final_val = computeIntegrandCore(v, h1, h2, t1, t2, cc_vec, dd_sqrt2pi_inv_vec, two_dd_sq_inv_vec, knot_x, omega, pp_omega, pZero, q0, qNA, C_pZero, C_q0, C_qNA, VAILD, ZERO, R_CENSORED, L_CENSORED, INVAILD, is_na1, is_na2, inv_SQRT2, inv_SQRT2PI)

    p_v = (v.^(-1) .* dd_sqrt2pi_inv_vec) .* exp(-(log(v) - cc_vec).^2 .* two_dd_sq_inv_vec);

    if v <= knot_x(1)
        inv_om_v = 1 / omega(1);
    elseif v >= knot_x(end)
        inv_om_v = 1 / omega(end);
    else
        inv_om_v = 1 / ppval(pp_omega, v);
    end

    inv_om_v_SQRT2PI = inv_om_v * inv_SQRT2PI;

    if (t1 == INVAILD) || (t2 == INVAILD)
        if t1 == INVAILD
            valid_h = h2; vt = t2; is_na = is_na2;
        else
            valid_h = h1; vt = t1; is_na = is_na1;
        end

        if is_na
            v_comp = -ones(size(valid_h));
        else
            z = (valid_h - v) .* inv_om_v;
            if vt == R_CENSORED
                v_comp = 0.5 .* erfc(z .* inv_SQRT2);
            elseif vt == L_CENSORED
                v_comp = 0.5 .* (1 + erf(z .* inv_SQRT2));
            else
                v_comp = exp(-0.5 .* (z.^2)) .* inv_om_v_SQRT2PI;
            end
        end

        if (vt == ZERO || vt == L_CENSORED)
            core = qNA .* (pZero + C_pZero .* (q0 + C_q0 .* v_comp));
        else
            core = qNA .* (C_pZero .* C_q0 .* v_comp);
        end
    else
        if is_na1
            v1 = -ones(size(h1));
        else
            z = (h1 - v) .* inv_om_v;
            if t1 == R_CENSORED
                v1 = 0.5 .* erfc(z .* inv_SQRT2);
            elseif t1 == L_CENSORED
                v1 = 0.5 .* (1 + erf(z .* inv_SQRT2));
            else
                v1 = exp(-0.5 .* (z.^2)) .* inv_om_v_SQRT2PI;
            end
        end

        if is_na2
            v2 = -ones(size(h2));
        else
            z = (h2 - v) .* inv_om_v;
            if t2 == R_CENSORED
                v2 = 0.5 .* erfc(z .* inv_SQRT2);
            elseif t2 == L_CENSORED
                v2 = 0.5 .* (1 + erf(z .* inv_SQRT2));
            else
                v2 = exp(-0.5 .* (z.^2)) .* inv_om_v_SQRT2PI;
            end
        end

        if (t1 == VAILD && t2 == VAILD)
            core = C_qNA .* C_pZero .* ((C_q0 .* v1) .* (C_q0 .* v2));
        elseif (t1 == ZERO && t2 == ZERO) || (t1 == L_CENSORED && t2 == L_CENSORED)
            core = C_qNA .* (pZero + C_pZero .* ((q0 + C_q0 .* v1) .* (q0 + C_q0 .* v2)));
        else
            if (t1==ZERO || t1==L_CENSORED), term1 = q0 + C_q0 .* v1; else, term1 = C_q0 .* v1; end
            if (t2==ZERO || t2==L_CENSORED), term2 = q0 + C_q0 .* v2; else, term2 = C_q0 .* v2; end
            core = C_qNA .* C_pZero .* (term1 .* term2);
        end
    end

    final_val = core .* p_v;
end
