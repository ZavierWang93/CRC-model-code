function parameters=transformParameters(paras,mc)

    last_index=0;

    prevalence=paras(:,1:(mc.num_prevalence-1));
    parameters.prevalence=getPrevalence(prevalence);
    last_index = last_index+size(parameters.prevalence,2)-1;

    parameters.ascertainment_ratio=paras(:,last_index+(1:mc.num_ascertainment_ratio));
    last_index = last_index+size(parameters.ascertainment_ratio,2);

    parameters.inc_NA=paras(:,last_index+(1:mc.num_inc_NA));
    last_index = last_index+size(parameters.inc_NA,2);

    parameters.inc_LSL=paras(:,last_index+(1:mc.num_inc_LSL));
    last_index = last_index+size(parameters.inc_LSL,2);

    parameters.exit_rate=paras(:,last_index+(1:mc.num_exit_rate));
    last_index = last_index+size(parameters.exit_rate,2);

    parameters.age_effect_on_disease_progression=paras(:,last_index+(1:mc.num_age_effect_on_disease_progression))*10;
    last_index = last_index+size(parameters.age_effect_on_disease_progression,2);

    parameters.prob_diagnosed=paras(:,last_index+(1:mc.num_prob_diagnosed));
    last_index = last_index+size(parameters.prob_diagnosed,2);

    parameters.mean_v{1}=paras(:,last_index+(1:mc.num_mean_v))*10;
    last_index = last_index+size(parameters.mean_v{1},2);

    parameters.std_v{1}=paras(:,last_index+(1:mc.num_std_v))*10;
    last_index = last_index+size(parameters.std_v{1},2);

    parameters.pZero{1}=paras(:,last_index+(1:mc.num_pZero));
    last_index = last_index+size(parameters.pZero{1},2);

    parameters.q0{1}=paras(:,last_index+(1:mc.num_q0));
    last_index = last_index+size(parameters.q0{1},2);

    parameters.qNA{1}=paras(:,last_index+(1:mc.num_qNA));
    last_index = last_index+size(parameters.qNA{1},2);

    parameters.omega{1}=paras(:,last_index+(1:mc.num_omega{1}))*10;
    last_index = last_index+size(parameters.omega{1},2);

    parameters.age_eff_mean{1}=paras(:,last_index+(1:mc.num_age_eff_mean))*10;
    last_index = last_index+size(parameters.age_eff_mean{1},2);

    parameters.mean_v{2}=paras(:,last_index+(1:mc.num_mean_v))*10;
    last_index = last_index+size(parameters.mean_v{2},2);

    parameters.std_v{2}=paras(:,last_index+(1:mc.num_std_v))*10;
    last_index = last_index+size(parameters.std_v{2},2);

    parameters.pZero{2}=paras(:,last_index+(1:mc.num_pZero));
    last_index = last_index+size(parameters.pZero{2},2);

    parameters.q0{2}=paras(:,last_index+(1:mc.num_q0));
    last_index = last_index+size(parameters.q0{2},2);

    parameters.qNA{2}=paras(:,last_index+(1:mc.num_qNA));
    last_index = last_index+size(parameters.qNA{2},2);

    parameters.omega{2}=paras(:,last_index+(1:mc.num_omega{2}))*10;
    last_index = last_index+size(parameters.omega{2},2);

    parameters.age_eff_mean{2}=paras(:,last_index+(1:mc.num_age_eff_mean))*10;
    last_index = last_index+size(parameters.age_eff_mean{2},2);

    parameters.knot_omega{1}=paras(:,last_index+(1:mc.num_knot_omega{1}));
    last_index = last_index+size(parameters.knot_omega{1},2);

    parameters.knot_omega{2}=paras(:,last_index+(1:mc.num_knot_omega{2}));
    last_index = last_index+size(parameters.knot_omega{2},2);

    parameters.mu=parameters.mean_v;
    parameters.sigma=parameters.std_v;

end
