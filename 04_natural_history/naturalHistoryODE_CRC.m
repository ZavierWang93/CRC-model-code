function dydt=naturalHistoryODE_CRC(tt,yy,Trans_A,all_cause_mortality,paras,mc)

    NORMAL=mc.NORMAL;
    NA=mc.NA;
    LSL=mc.LSL;
    AA=mc.AA;
    HSL=mc.HSL;
    CRC_UNDIAGNOSED=mc.CRC_UNDIAGNOSED;
    CRC_SL=mc.CRC_SL;
    CRC_SURVIVOR=mc.CRC_SURVIVOR;
    DEATH=mc.DEATH;
    age_effect_progression=[1 paras.age_effect_on_disease_progression];

    if (tt<mc.knot_age_effect_disease_progression(1))
        age_effect=1;
    else
        age_effect=pchip(mc.knot_age_effect_disease_progression,age_effect_progression,tt);
    end

    age_effect=max(0,age_effect);

    for stage=[AA,HSL,CRC_UNDIAGNOSED(1),CRC_UNDIAGNOSED(2),CRC_UNDIAGNOSED(3),CRC_UNDIAGNOSED(4)]
        Trans_A(:, stage) = Trans_A(:, stage)*age_effect;
    end

    if tt>mc.knot_incidence_rate(end)
        incidence_rate_AA=paras.inc_NA(end);
        incidence_rate_LSL=paras.inc_LSL(end);

    else
        incidence_rate_AA=pchip(mc.knot_incidence_rate,paras.inc_NA,tt);
        incidence_rate_LSL=pchip(mc.knot_incidence_rate,paras.inc_LSL,tt);

    end

    num_stages = numel(yy);

    Trans_B = zeros(num_stages);

    Trans_B(NA,NORMAL)=incidence_rate_AA;

    Trans_B(LSL,NORMAL)=incidence_rate_LSL;

    Trans_B(NORMAL,NORMAL)=-incidence_rate_AA-incidence_rate_LSL;

    death_rate = interp1(mc.age_knots_HKCR,all_cause_mortality,tt);

    for stage=NORMAL:CRC_SURVIVOR
        Trans_B(stage, stage) = Trans_B(stage, stage)-death_rate;
        Trans_B(DEATH, stage) = Trans_B(DEATH, stage)+death_rate;
    end

    dydt=(Trans_A+Trans_B)*yy;

end