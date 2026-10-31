function output=logLikelihood_CRC(initial_para_M,initial_para_F,observe_sorts_M,observe_sorts_F,mc,inm_M,inm_F)

    safe_log = @safeLog;

    use_prev_cache = isfield(mc, 'copy_use_prev_cache') && mc.copy_use_prev_cache;
    prev_cache = [];
    if use_prev_cache && isfield(mc, 'copy_prev_cache')
        prev_cache = mc.copy_prev_cache;
    end

    if ~isempty(prev_cache) && (~isfield(prev_cache, 'state_names') || ...
            ~strcmp(prev_cache.state_names, 'LSL_HSL_pathSL') || ~isfield(prev_cache, 'crc_stage_reporting') || ...
            ~strcmp(prev_cache.crc_stage_reporting, 'asrstage_v1'))
        prev_cache = [];
    end

    parameter_count = size(initial_para_M, 2);
    parameters = syncMCMCParameters([initial_para_M, initial_para_F], parameter_count);
    initial_para_M = parameters(:, 1:parameter_count);
    initial_para_F = parameters(:, parameter_count+1:end);
    paras_cell = {initial_para_M, initial_para_F};
    inm_cell = {inm_M, inm_F};
    obs_cell = {observe_sorts_M, observe_sorts_F};

    MODEL = mc.MODEL;
    OBSERVED = mc.OBSERVED;
    CRC_INCIDENCE = mc.CRC_INCIDENCE;
    CRC_MORTALITY = mc.CRC_MORTALITY;
    NORMAL = mc.NORMAL;
    UNKNOWN = mc.UNKNOWN;
    CRC = mc.CRC;
    CRC_UNDIAGNOSED = mc.CRC_UNDIAGNOSED;
    col_sens = mc.colonoscopy_sensitivity;
    No_CS_pathway_pro=mc.No_CS_pathway_pro;

    cohort_cell = cell(1, 2);
    logL_base = zeros(1, 2);
    prob_disease_stage_cell = cell(1, 2);
    model_specs_cell = cell(2, max(mc.contractor_list));

    same_para_gender = false(1, 2);
    ode_changed_gender = true(1, 2);
    ode_likelihood_changed_gender = true(1, 2);
    fit_changed_gender = true(1, 2);
    reuse_ode_gender = false(1, 2);
    reuse_task_gender = false(1, 2);
    fit_changed_gender_lab = true(2, max(mc.contractor_list));

    num_para = size(paras_cell{1}, 2);

    n_prev = mc.num_prevalence - 1;
    n_asr = mc.num_ascertainment_ratio;
    n_inc_na = mc.num_inc_NA;
    n_inc_lsl = mc.num_inc_LSL;
    n_exit = mc.num_exit_rate;
    n_ageeff = mc.num_age_effect_on_disease_progression;
    n_diag = mc.num_prob_diagnosed;

    idx_prev = 1:n_prev;
    idx_asr = (n_prev + 1):(n_prev + n_asr);
    k = n_prev + n_asr;
    idx_inc_na = (k + 1):(k + n_inc_na); k = k + n_inc_na;
    idx_inc_lsl = (k + 1):(k + n_inc_lsl); k = k + n_inc_lsl;
    idx_exit = (k + 1):(k + n_exit); k = k + n_exit;
    idx_ageeff = (k + 1):(k + n_ageeff); k = k + n_ageeff;
    idx_diag = (k + 1):(k + n_diag);

    default_ode_idx = [idx_prev, idx_asr, idx_inc_na, idx_inc_lsl, idx_exit, idx_ageeff, idx_diag];
    default_ode_likelihood_idx = [idx_prev, idx_inc_na, idx_inc_lsl, idx_exit, idx_ageeff, idx_diag];

    if isfield(mc, 'copy_ode_param_idx') && ~isempty(mc.copy_ode_param_idx)
        ode_param_idx = unique(mc.copy_ode_param_idx(:)');
    else
        ode_param_idx = default_ode_idx;
    end

    if isfield(mc, 'copy_ode_likelihood_param_idx') && ~isempty(mc.copy_ode_likelihood_param_idx)
        ode_likelihood_param_idx = unique(mc.copy_ode_likelihood_param_idx(:)');
    else
        ode_likelihood_param_idx = default_ode_likelihood_idx;
    end

    if isfield(mc, 'copy_fit_param_idx') && ~isempty(mc.copy_fit_param_idx)
        fit_param_idx = unique(mc.copy_fit_param_idx(:)');
    else
        fit_param_idx = setdiff(1:num_para, ode_param_idx);
    end

    ode_param_idx = ode_param_idx(ode_param_idx >= 1 & ode_param_idx <= num_para);
    ode_likelihood_param_idx = ode_likelihood_param_idx(ode_likelihood_param_idx >= 1 & ode_likelihood_param_idx <= num_para);
    fit_param_idx = fit_param_idx(fit_param_idx >= 1 & fit_param_idx <= num_para);

    ode_param_idx = unique([ode_param_idx, idx_asr]);
    ode_likelihood_param_idx = unique([ode_likelihood_param_idx, idx_asr]);

    scenario_config = struct('kind', 'nh_meanage_base', 'fit_abs_tol', 1e-10, 'fit_rel_tol', 1e-6);
    cache_scenario_matches = ~isempty(prev_cache) && ...
        isfield(prev_cache, 'scenario_config') && ...
        isequaln(prev_cache.scenario_config, scenario_config);
    if use_prev_cache && cache_scenario_matches && isfield(prev_cache, 'paras_cell') && ...
            iscell(prev_cache.paras_cell) && numel(prev_cache.paras_cell) == 2 && ...
            isfield(prev_cache, 'cohort_cell') && iscell(prev_cache.cohort_cell) && ...
            numel(prev_cache.cohort_cell) == 2 && ...
            isfield(prev_cache, 'logL_base') && numel(prev_cache.logL_base) == 2 && ...
            isfield(prev_cache, 'prob_disease_stage_cell') && ...
            iscell(prev_cache.prob_disease_stage_cell) && ...
            numel(prev_cache.prob_disease_stage_cell) == 2
        for g = 1:2
            par_prev = prev_cache.paras_cell{g};
            par_curr = paras_cell{g};
            same_para_gender(g) = isequal(par_prev, par_curr);

            if ~isempty(ode_param_idx)
                ode_changed_gender(g) = ~isequal(par_prev(:, ode_param_idx), par_curr(:, ode_param_idx));
            else
                ode_changed_gender(g) = false;
            end

            if ~isempty(ode_likelihood_param_idx)
                ode_likelihood_changed_gender(g) = ~isequal(par_prev(:, ode_likelihood_param_idx), par_curr(:, ode_likelihood_param_idx));
            else
                ode_likelihood_changed_gender(g) = false;
            end

            if ~isempty(fit_param_idx)
                fit_changed_gender(g) = ~isequal(par_prev(:, fit_param_idx), par_curr(:, fit_param_idx));
            else
                fit_changed_gender(g) = false;
            end

            reuse_ode_gender(g) = ~ode_changed_gender(g);
            reuse_task_gender(g) = ~(ode_changed_gender(g) || fit_changed_gender(g));
        end
    end

    for g = 1:2
        if reuse_ode_gender(g)
            cohort_cell{g} = prev_cache.cohort_cell{g};
            logL_base(g) = prev_cache.logL_base(g);
            prob_disease_stage_cell{g} = prev_cache.prob_disease_stage_cell{g};
        else
            cohort_cell{g} = simulateNaturalHistory_CRC(g, paras_cell{g}, mc);
            cohort_cell{g}.stats_5y{CRC_INCIDENCE, OBSERVED} = inm_cell{g}.stats_5y{CRC_INCIDENCE, OBSERVED};
            cohort_cell{g}.stats_5y{CRC_MORTALITY, OBSERVED} = inm_cell{g}.stats_5y{CRC_MORTALITY, OBSERVED};

            for jj = 1:5
                obs_val = cohort_cell{g}.stats_5y{CRC_INCIDENCE, OBSERVED}(:,jj);
                exp_val = max(1e-300, cohort_cell{g}.stats_5y{CRC_INCIDENCE, MODEL}(:,jj)/No_CS_pathway_pro);
                logL_base(g) = logL_base(g) + obs_val'*log(exp_val) - sum(exp_val);
            end

            obs_val = cohort_cell{g}.stats_5y{CRC_MORTALITY, OBSERVED};
            exp_val = max(1e-300, cohort_cell{g}.stats_5y{CRC_MORTALITY, MODEL}/No_CS_pathway_pro);
            logL_base(g) = logL_base(g) + obs_val*log(exp_val') - sum(exp_val);

            prob_disease_stage_cell{g} = cohort_cell{g}.prevalence_1y(:, NORMAL:CRC_UNDIAGNOSED(4), 2);

        end

        for lab = mc.contractor_list
            model_specs_cell{g, lab} = transformModelSpecs(paras_cell{g}, mc, lab);
        end
    end

    if use_prev_cache && ~isempty(prev_cache) && isfield(prev_cache, 'model_specs_cell')
        for g = 1:2
            for lab = mc.contractor_list
                ms_prev_ok = size(prev_cache.model_specs_cell,1) >= g && size(prev_cache.model_specs_cell,2) >= lab;
                if ms_prev_ok && ~isempty(prev_cache.model_specs_cell{g, lab})
                    fit_changed_gender_lab(g, lab) = ~isequal(prev_cache.model_specs_cell{g, lab}, model_specs_cell{g, lab});
                else
                    fit_changed_gender_lab(g, lab) = true;
                end
            end
        end

        for g = 1:2
            fit_changed_gender(g) = any(fit_changed_gender_lab(g, mc.contractor_list));
            reuse_task_gender(g) = ~(ode_likelihood_changed_gender(g) || fit_changed_gender(g));
        end
    else
        for g = 1:2
            for lab = mc.contractor_list
                fit_changed_gender_lab(g, lab) = fit_changed_gender(g);
            end
        end
    end

    task_gender = []; task_lab = []; task_outcomes = []; task_stages = [];
    task_data = {}; task_ms = {};

    for g = 1:2
        for lab = mc.contractor_list
            for outcomes = NORMAL:UNKNOWN
                if (outcomes == UNKNOWN || outcomes == NORMAL)
                    possible_stages = NORMAL:CRC;
                else
                    possible_stages = outcomes;
                end
                for stages = possible_stages

                    data_candidate = obs_cell{g}{outcomes, lab};
                    has_data = false;
                    for ci = 1:numel(data_candidate)
                        if ~isempty(data_candidate{ci})
                            has_data = true;
                            break;
                        end
                    end
                    if ~has_data
                        continue;
                    end

                    ms_full = model_specs_cell{g, lab};
                    ms_slim.knot_omega = ms_full.knot_omega;
                    ms_slim.omega = ms_full.omega;
                    ms_slim.pZero = ms_full.pZero;
                    ms_slim.qNA = ms_full.qNA;
                    ms_slim.q0 = ms_full.q0;
                    ms_slim.mu = ms_full.mu;
                    ms_slim.sigma = ms_full.sigma;
                    ms_slim.age_eff_mean = ms_full.age_eff_mean;

                    task_gender(end+1) = g;
                    task_lab(end+1) = lab;
                    task_outcomes(end+1) = outcomes;
                    task_stages(end+1) = stages;
                    task_data{end+1} = data_candidate;
                    task_ms{end+1} = ms_slim;
                end
            end
        end
    end

    num_tasks = length(task_gender);

    par_prob = cell(num_tasks, 1);
    par_age = cell(num_tasks, 1);
    par_id = cell(num_tasks, 1);
    par_crc_stage = cell(num_tasks, 1);
    par_likelihood_col = cell(num_tasks, 1);
    par_profile = cell(num_tasks, 1);

    reuse_task = false(num_tasks, 1);
    prev_task_cache_ok = false;
    if use_prev_cache && ~isempty(prev_cache) && isfield(prev_cache, 'task_gender') && ...
            isfield(prev_cache, 'task_lab') && isfield(prev_cache, 'task_outcomes') && ...
            isfield(prev_cache, 'task_stages') && ...
            isfield(prev_cache, 'par_prob') && isfield(prev_cache, 'par_age') && ...
            isfield(prev_cache, 'par_id') && isfield(prev_cache, 'par_crc_stage') && ...
            isfield(prev_cache, 'par_likelihood_col') && isfield(prev_cache, 'par_numpts') && ...
            isfield(prev_cache, 'par_profile')
        prev_task_cache_ok = numel(prev_cache.task_gender) == num_tasks && ...
            isequal(prev_cache.task_gender(:), task_gender(:)) && ...
            isequal(prev_cache.task_lab(:), task_lab(:)) && ...
            isequal(prev_cache.task_outcomes(:), task_outcomes(:)) && ...
            isequal(prev_cache.task_stages(:), task_stages(:));
        if prev_task_cache_ok
            for task_idx = 1:num_tasks
                g = task_gender(task_idx);
                lab = task_lab(task_idx);
                reuse_task(task_idx) = (~ode_likelihood_changed_gender(g)) && (~fit_changed_gender_lab(g, lab));
            end
        end
    end

    age_1y_M = cohort_cell{1}.age_1y;
    age_1y_F = cohort_cell{2}.age_1y;
    pds_M = prob_disease_stage_cell{1};
    pds_F = prob_disease_stage_cell{2};
    crc_undiag_vec = CRC_UNDIAGNOSED;
    asr_M = reshape(cohort_cell{1}.ascertainment_ratio_this, [], numel(crc_undiag_vec));
    asr_F = reshape(cohort_cell{2}.ascertainment_ratio_this, [], numel(crc_undiag_vec));

    asr_M(end+1, :) = paras_cell{1}(idx_asr([3 6 9 9]));
    asr_F(end+1, :) = paras_cell{2}(idx_asr([3 6 9 9]));

    mc_slim_base.ZERO = mc.ZERO;
    mc_slim_base.EXACT = mc.EXACT;
    mc_slim_base.INVALID = mc.INVALID;
    mc_slim_base.INVALID_FIT = mc.INVALID_FIT;
    mc_slim_base.RIGHT_CENSORED = mc.RIGHT_CENSORED;
    mc_slim_base.LEFT_CENSORED = mc.LEFT_CENSORED;
    mc_slim_base.logFIT_left_censored = mc.logFIT_left_censored;
    mc_slim_base.logFIT_right_censored = mc.logFIT_right_censored;
    mc_slim_base.CRC = mc.CRC;
    if isfield(mc, 'dyn_group_profile')
        mc_slim_base.dyn_group_profile = mc.dyn_group_profile;
    end

    has_mem = isfield(mc, 'group_numpts_start') && iscell(mc.group_numpts_start);

    par_numpts = cell(num_tasks, 1);
    parfor task_idx = 1:num_tasks
        g = task_gender(task_idx);
        lab = task_lab(task_idx);
        outcomes = task_outcomes(task_idx);
        stages = task_stages(task_idx);

        if reuse_task(task_idx) && prev_task_cache_ok
            par_prob{task_idx} = prev_cache.par_prob{task_idx};
            par_age{task_idx} = prev_cache.par_age{task_idx};
            par_id{task_idx} = prev_cache.par_id{task_idx};
            par_crc_stage{task_idx} = prev_cache.par_crc_stage{task_idx};
            par_likelihood_col{task_idx} = prev_cache.par_likelihood_col{task_idx};
            par_numpts{task_idx} = prev_cache.par_numpts{task_idx};
            par_profile{task_idx} = prev_cache.par_profile{task_idx};
            continue;
        end

        data_this = task_data{task_idx};
        ms_this = task_ms{task_idx};

        mc_task = mc_slim_base;
        if has_mem
            if g <= size(mc.group_numpts_start, 1) && ...
               lab <= size(mc.group_numpts_start, 2) && ...
               outcomes <= size(mc.group_numpts_start, 3) && ...
               stages <= size(mc.group_numpts_start, 4)
                mc_task.group_numpts_start = mc.group_numpts_start{g, lab, outcomes, stages};
            end
        end

        output_this = integrateFITLikelihood( ...
            data_this, stages, ms_this, mc_task, lab);

        if g == 1
            age_1y = age_1y_M; pds = pds_M;
            asr = asr_M;
        else
            age_1y = age_1y_F; pds = pds_F;
            asr = asr_F;
        end

        prevalence_this = zeros(length(output_this.age), 1);
        if outcomes == CRC
            for stage_idx = crc_undiag_vec
                fI = find(output_this.CRC_stage == stage_idx - crc_undiag_vec(1) + 1);
                if ~isempty(fI)

                    stage_num = stage_idx - crc_undiag_vec(1) + 1;
                    prevalence_this(fI, 1) = interp1(age_1y, pds(:, stage_idx), output_this.age(fI)) .* ...
                        interp1(age_1y, asr(:, stage_num), output_this.age(fI));
                end
            end
            fI = find(output_this.CRC_stage == 5);
            if ~isempty(fI)

                stage_prevalence = interp1(age_1y, pds(:, crc_undiag_vec), output_this.age(fI));
                stage_asr = interp1(age_1y, asr, output_this.age(fI));
                prevalence_this(fI, 1) = sum(stage_prevalence .* (1 - stage_asr), 2);
            end
        else
            if stages == CRC
                prevalence_this = interp1(age_1y, sum(pds(:, crc_undiag_vec), 2), output_this.age);
            else
                prevalence_this = interp1(age_1y, pds(:, stages), output_this.age);
            end
        end

        prob_FIT_this = output_this.prob;

        if (outcomes == NORMAL)
            if (stages == NORMAL)
                L_col = prevalence_this .* prob_FIT_this;
            else
                L_col = (1 - col_sens(stages)) * prevalence_this .* prob_FIT_this;
            end
        elseif (outcomes == UNKNOWN)
            L_col = prevalence_this .* prob_FIT_this;
        else
            L_col = col_sens(stages) * prevalence_this .* prob_FIT_this;
        end

        par_prob{task_idx} = output_this.prob;
        par_age{task_idx} = output_this.age;
        par_id{task_idx} = output_this.screenee_id;
        if isfield(output_this, 'CRC_stage')
            par_crc_stage{task_idx} = output_this.CRC_stage;
        end
        par_likelihood_col{task_idx} = L_col;
        if isfield(output_this, 'final_numpts_mat')
            par_numpts{task_idx} = output_this.final_numpts_mat;
        end
        if isfield(output_this, 'profile')
            par_profile{task_idx} = output_this.profile;
        end
    end

    like_ind_cell = {cell(UNKNOWN, max(mc.contractor_list)), cell(UNKNOWN, max(mc.contractor_list))};
    logL_comp_cell = {cell(UNKNOWN, max(mc.contractor_list)), cell(UNKNOWN, max(mc.contractor_list))};
    prob_log_FIT_obx = cell(1, 2);

    for task_idx = 1:num_tasks
        g = task_gender(task_idx);
        lab = task_lab(task_idx);
        outcomes = task_outcomes(task_idx);
        stages = task_stages(task_idx);

        out_slim.prob = par_prob{task_idx};
        out_slim.age = par_age{task_idx};
        out_slim.screenee_id = par_id{task_idx};
        if ~isempty(par_crc_stage{task_idx})
            out_slim.CRC_stage = par_crc_stage{task_idx};
        end

        prob_log_FIT_obx{g}{outcomes, stages, lab} = out_slim;

        needed_rows = length(par_likelihood_col{task_idx});
        current_rows = size(like_ind_cell{g}{outcomes, lab}, 1);
        if current_rows < needed_rows && current_rows > 0
            like_ind_cell{g}{outcomes, lab}(current_rows+1:needed_rows, :) = 0;
        end
        like_ind_cell{g}{outcomes, lab}(1:needed_rows, stages) = par_likelihood_col{task_idx};
    end

    logL_FIT = zeros(1, 2);
    for g = 1:2
        for lab = mc.contractor_list
            for outcomes = NORMAL:UNKNOWN
                if isempty(like_ind_cell{g}{outcomes, lab})
                    continue;
                end
                if (outcomes == UNKNOWN || outcomes == NORMAL)
                    logL_comp_cell{g}{outcomes, lab} = safe_log(sum(like_ind_cell{g}{outcomes, lab}, 2));
                else
                    logL_comp_cell{g}{outcomes, lab}(:, outcomes) = safe_log(like_ind_cell{g}{outcomes, lab}(:, outcomes));
                end
                logL_FIT(g) = logL_FIT(g) + sum(logL_comp_cell{g}{outcomes, lab}, 'all');
            end
        end
    end

    LogL_total = 0;
    ans_out = cell(1, 2);

    for g = 1:2
        ms_final = transformParameters(paras_cell{g}, mc);

        ans_out{g}.para = paras_cell{g};
        ans_out{g}.model_specs = ms_final;
        ans_out{g}.cohort = cohort_cell{g};
        ans_out{g}.prob_log_FIT_obx = prob_log_FIT_obx{g};

        LogL_g = logL_base(g) + logL_FIT(g);
        LogL_g = LogL_g + ascertainmentRatioConstraint(ms_final, mc);
        LogL_g = LogL_g + logPriorSL(cohort_cell{g},mc);
        LogL_g = LogL_g + cdfConstraint(ms_final, mc);
        LogL_g = LogL_g + prevalenceConstraint(ms_final, mc);
        LogL_g = LogL_g + incidenceConstraint(ms_final, mc);
        LogL_g = LogL_g + ageProgressionConstraint(ms_final);
        LogL_g = LogL_g + logPriorNA(ms_final);

        ans_out{g}.LogL = LogL_g;
        LogL_total = LogL_total + LogL_g;
    end
    LogL_total = LogL_total + logPriorOmega(ms_final);

    final_res_cells = cell(2, max(mc.contractor_list), 7, 6);
    if has_mem
        final_res_cells = mc.group_numpts_start;
    end

    for task_idx = 1:num_tasks
        if ~isempty(par_numpts{task_idx})
            g = task_gender(task_idx);
            lab = task_lab(task_idx);
            outcomes = task_outcomes(task_idx);
            stages = task_stages(task_idx);

            if isempty(final_res_cells{g, lab, outcomes, stages})
                final_res_cells{g, lab, outcomes, stages} = par_numpts{task_idx};
            else
                final_res_cells{g, lab, outcomes, stages} = max(final_res_cells{g, lab, outcomes, stages}, par_numpts{task_idx});
            end
        end
    end

    output.group_numpts_start = final_res_cells;
    output.final_numpts_mat = final_res_cells;

    prof_mat = [];
    for task_idx = 1:num_tasks
        p = par_profile{task_idx};
        if isempty(p) || ~isfield(p, 'while_time')
            continue;
        end
        n_prof = numel(p.while_time);
        if n_prof == 0
            continue;
        end
        gcol = repmat(task_gender(task_idx), n_prof, 1);
        labcol = repmat(task_lab(task_idx), n_prof, 1);
        outcol = repmat(task_outcomes(task_idx), n_prof, 1);
        stagecol = repmat(task_stages(task_idx), n_prof, 1);

        t1 = p.target_1(:);
        t2 = p.target_2(:);
        nrows = p.num_rows(:);
        wtime = p.while_time(:);

        if isfield(p, 'num_iter')
            niter = p.num_iter(:);
        else
            niter = nan(n_prof,1);
        end
        if isfield(p, 'final_numpts')
            fnum = p.final_numpts(:);
        else
            fnum = nan(n_prof,1);
        end

        prof_mat = [prof_mat; [gcol, labcol, outcol, stagecol, t1, t2, nrows, wtime, niter, fnum]];
    end
    output.profile_header = {'gender','lab','outcome','stage','target_1','target_2','num_rows','while_time','num_iter','final_numpts'};
    output.profile_matrix = prof_mat;

    cache_out.state_names = 'LSL_HSL_pathSL';
    cache_out.scenario_config = scenario_config;
    cache_out.crc_stage_reporting = 'asrstage_v1';
    cache_out.paras_cell = paras_cell;
    cache_out.cohort_cell = cohort_cell;
    cache_out.logL_base = logL_base;
    cache_out.prob_disease_stage_cell = prob_disease_stage_cell;
    cache_out.model_specs_cell = model_specs_cell;
    cache_out.task_gender = task_gender;
    cache_out.task_lab = task_lab;
    cache_out.task_outcomes = task_outcomes;
    cache_out.task_stages = task_stages;
    cache_out.par_prob = par_prob;
    cache_out.par_age = par_age;
    cache_out.par_id = par_id;
    cache_out.par_crc_stage = par_crc_stage;
    cache_out.par_likelihood_col = par_likelihood_col;
    cache_out.par_numpts = par_numpts;
    cache_out.par_profile = par_profile;

    output.copy_prev_cache = cache_out;
    output.copy_reuse_info.ode_reused_gender = reuse_ode_gender;
    output.copy_reuse_info.task_reused_gender = reuse_task_gender;
    output.copy_reuse_info.ode_changed_gender = ode_changed_gender;
    output.copy_reuse_info.ode_likelihood_changed_gender = ode_likelihood_changed_gender;
    output.copy_reuse_info.fit_changed_gender = fit_changed_gender;
    output.copy_reuse_info.fit_changed_gender_lab = fit_changed_gender_lab;
    output.copy_reuse_info.same_para_gender = same_para_gender;
    output.copy_reuse_info.ode_param_idx = ode_param_idx;
    output.copy_reuse_info.ode_likelihood_param_idx = ode_likelihood_param_idx;
    output.copy_reuse_info.fit_param_idx = fit_param_idx;
    output.copy_reuse_info.reused_tasks = sum(reuse_task);
    output.copy_reuse_info.total_tasks = num_tasks;

    output.LogL = LogL_total;
    output.ans{1} = ans_out{1};
    output.ans{2} = ans_out{2};
end
