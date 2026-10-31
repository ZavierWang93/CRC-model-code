function model_specs_use=transformModelSpecs(paras,mc,lab)

    model_specs=transformParameters(paras,mc);

    model_specs_use.age_effect_on_disease_progression=model_specs.age_effect_on_disease_progression;

    model_specs_use.ascertainment_ratio=model_specs.ascertainment_ratio;

    model_specs_use.exit_rate=model_specs.exit_rate;

    model_specs_use.inc_NA=model_specs.inc_NA;

    model_specs_use.inc_LSL=model_specs.inc_LSL;

    model_specs_use.knot_omega=model_specs.knot_omega{lab};

    model_specs_use.mu=model_specs.mean_v{lab};

    model_specs_use.omega=model_specs.omega{lab};

    model_specs_use.prevalence=model_specs.prevalence;

    model_specs_use.prob_diagnosed=model_specs.prob_diagnosed;

    model_specs_use.pZero=model_specs.pZero{lab};

    model_specs_use.q0=model_specs.q0{lab};

    model_specs_use.qNA=model_specs.qNA{lab};

    model_specs_use.sigma=model_specs.std_v{lab};

    model_specs_use.age_eff_mean=model_specs.age_eff_mean{lab};

end