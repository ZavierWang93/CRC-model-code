function cohort=simulateNaturalHistory_CRC(gender_this,paras,mc)

parameters=transformParameters(paras,mc);
cohort.parameters=paras;
prevalence_ini=parameters.prevalence/sum(parameters.prevalence);

cohort_rates=getRates(gender_this,mc);
all_cause_mortality=cohort_rates(gender_this).interpolated_rate{mc.ALL_CAUSE_MORTALITY};

Initial_population_1=mc.cohort_size.*prevalence_ini;
Initial_population=[Initial_population_1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

NORMAL=mc.NORMAL;
NA=mc.NA;
LSL=mc.LSL;
AA=mc.AA;
HSL=mc.HSL;
CRC=mc.CRC;
CRC_UNDIAGNOSED=mc.CRC_UNDIAGNOSED;
CRC_DIAGNOSED=mc.CRC_DIAGNOSED;
CRC_SURVIVOR=mc.CRC_SURVIVOR;
CRC_DEATH=mc.CRC_DEATH;

CUM_CRC_DIAGNOSED=mc.CUM_CRC_DIAGNOSED;
CRC_AA=mc.CRC_AA;
CRC_SL=mc.CRC_SL;

age_range=mc.age_range;
age_knots_HKCR=mc.age_knots_HKCR;
MODEL=mc.MODEL;
OBSERVED=mc.OBSERVED;

COHORT_TOTAL=mc.COHORT_TOTAL;
COHORT_WITHOUT_DIAGNOSED_CRC=mc.COHORT_WITHOUT_DIAGNOSED_CRC;

CRC_INCIDENCE=mc.CRC_INCIDENCE;
CRC_MORTALITY=mc.CRC_MORTALITY;
CRC_INCIDENCE_RATE=mc.CRC_INCIDENCE_RATE;
CRC_MORTALITY_RATE=mc.CRC_MORTALITY_RATE;
ALL_CAUSE_MORTALITY=mc.ALL_CAUSE_MORTALITY;
ALL_CAUSE_MORTALITY_RATE=mc.ALL_CAUSE_MORTALITY_RATE;

exit_rate=parameters.exit_rate;
diagnosed_rate=parameters.prob_diagnosed;
CRCdeath=(1-mc.CRC_5y_survival_prob);

ascertainment_ratio_original=parameters.ascertainment_ratio;
ascertainment_ratio=[ascertainment_ratio_original ascertainment_ratio_original(7:9)];

num_CRC_stages=mc.num_CRC_stages;

Trans_A=zeros(mc.CRC_SL,mc.CRC_SL);

Trans_A(NA,NA)=-exit_rate(1);
Trans_A(AA,NA)=exit_rate(1);

Trans_A(LSL,LSL)=-exit_rate(3);
Trans_A(HSL,LSL)=exit_rate(3);

Trans_A(AA,AA)=-exit_rate(2);
Trans_A(CRC_UNDIAGNOSED(1),AA)=exit_rate(2);
Trans_A(CRC_AA,AA)=exit_rate(2);

Trans_A(HSL,HSL)=-(exit_rate(2)+exit_rate(4));
Trans_A(CRC_UNDIAGNOSED(1),HSL)=(exit_rate(2)+exit_rate(4));
Trans_A(CRC_SL,HSL)=(exit_rate(2)+exit_rate(4));

Trans_A(CRC_UNDIAGNOSED(1),CRC_UNDIAGNOSED(1))=-exit_rate(5);
Trans_A(CRC_DIAGNOSED(1),CRC_UNDIAGNOSED(1))=exit_rate(5)*diagnosed_rate(1);
Trans_A(CUM_CRC_DIAGNOSED(1),CRC_UNDIAGNOSED(1))=exit_rate(5)*diagnosed_rate(1);
Trans_A(CRC_UNDIAGNOSED(2),CRC_UNDIAGNOSED(1))=exit_rate(5)*(1-diagnosed_rate(1));

Trans_A(CRC_UNDIAGNOSED(2),CRC_UNDIAGNOSED(2))=-exit_rate(6);
Trans_A(CRC_DIAGNOSED(2),CRC_UNDIAGNOSED(2))=exit_rate(6)*diagnosed_rate(2);
Trans_A(CUM_CRC_DIAGNOSED(2),CRC_UNDIAGNOSED(2))=exit_rate(6)*diagnosed_rate(2);
Trans_A(CRC_UNDIAGNOSED(3),CRC_UNDIAGNOSED(2))=exit_rate(6)*(1-diagnosed_rate(2));

Trans_A(CRC_UNDIAGNOSED(3),CRC_UNDIAGNOSED(3))=-exit_rate(7);
Trans_A(CRC_DIAGNOSED(3),CRC_UNDIAGNOSED(3))=exit_rate(7)*diagnosed_rate(3);
Trans_A(CUM_CRC_DIAGNOSED(3),CRC_UNDIAGNOSED(3))=exit_rate(7)*diagnosed_rate(3);
Trans_A(CRC_UNDIAGNOSED(4),CRC_UNDIAGNOSED(3))=exit_rate(7)*(1-diagnosed_rate(3));

Trans_A(CRC_UNDIAGNOSED(4),CRC_UNDIAGNOSED(4))=-exit_rate(8);
Trans_A(CRC_DIAGNOSED(4),CRC_UNDIAGNOSED(4))=exit_rate(8)*1;
Trans_A(CUM_CRC_DIAGNOSED(4),CRC_UNDIAGNOSED(4))=exit_rate(8)*1;

rate=-log(0.001)/5;
Trans_A(CRC_DIAGNOSED(1),CRC_DIAGNOSED(1))=0-rate;
Trans_A(CRC_DEATH,CRC_DIAGNOSED(1))=rate*CRCdeath(1);
Trans_A(CRC_SURVIVOR,CRC_DIAGNOSED(1))=rate-rate*CRCdeath(1);

Trans_A(CRC_DIAGNOSED(2),CRC_DIAGNOSED(2))=0-rate;
Trans_A(CRC_DEATH,CRC_DIAGNOSED(2))=rate*CRCdeath(2);
Trans_A(CRC_SURVIVOR,CRC_DIAGNOSED(2))=rate-rate*CRCdeath(2);

Trans_A(CRC_DIAGNOSED(3),CRC_DIAGNOSED(3))=0-rate;
Trans_A(CRC_DEATH,CRC_DIAGNOSED(3))=rate*CRCdeath(3);
Trans_A(CRC_SURVIVOR,CRC_DIAGNOSED(3))=rate-rate*CRCdeath(3);

Trans_A(CRC_DIAGNOSED(4),CRC_DIAGNOSED(4))=0-rate;
Trans_A(CRC_DEATH,CRC_DIAGNOSED(4))=rate*CRCdeath(4);
Trans_A(CRC_SURVIVOR,CRC_DIAGNOSED(4))=rate-rate*CRCdeath(4);

[sln_time, sln_state] = ode45(@(tt,yy)naturalHistoryODE_CRC(tt, yy, Trans_A, all_cause_mortality, parameters, mc), ...
    age_range, Initial_population);

cohort.state_sln = sln_state;

cohort.age_sln = sln_time;
cohort.size_sln = zeros(numel(sln_time),2);
cohort.prevalence_sln = zeros(numel(sln_time), CRC_SURVIVOR, 2);

for idx = [COHORT_TOTAL COHORT_WITHOUT_DIAGNOSED_CRC]
        if (idx==1)
            states_included = NORMAL:CRC_SURVIVOR;
        else
            states_included = NORMAL:CRC_UNDIAGNOSED(end);
        end
        cohort.size_sln(:, idx) = sum(cohort.state_sln(:, states_included), 2);
        cohort.prevalence_sln(:, states_included, idx) = cohort.state_sln(:, states_included)./repmat(cohort.size_sln(:, idx), 1, numel(states_included));
end

cohort.age_1y = age_range(1):age_range(end);
cohort.state_1y = interp1(sln_time, sln_state,cohort.age_1y);
cohort.size_1y = interp1(cohort.age_sln, cohort.size_sln, cohort.age_1y);
cohort.prevalence_1y = interp1(cohort.age_sln, cohort.prevalence_sln, cohort.age_1y);

    xx = cohort.state_1y(:, CRC_UNDIAGNOSED(1:end));
    cohort.undiagnosed_CRC_stage_dist_1y = xx./repmat(sum(xx, 2), 1, size(xx, 2));

    for CRC_stage = 1:mc.num_CRC_stages
        cohort.incidence_CRC_1y(:, CRC_stage) = diff(cohort.state_1y(:, CUM_CRC_DIAGNOSED(CRC_stage)));
    end

    cohort.incidence_CRC_SL_1y = diff(cohort.state_1y(:, CRC_SL));
    cohort.incidence_CRC_AA_1y = diff(cohort.state_1y(:, CRC_AA));
    cohort.age_specific_prop_CRC_SL = cohort.incidence_CRC_SL_1y./(cohort.incidence_CRC_SL_1y+cohort.incidence_CRC_AA_1y);
    cohort.overall_prop_CRC_SL = sum(cohort.incidence_CRC_SL_1y)./sum(cohort.incidence_CRC_SL_1y+cohort.incidence_CRC_AA_1y);

age_vec = age_range(1):age_range(end);
age_vec = age_vec(1:size(cohort.incidence_CRC_1y, 1));

cohort.incidence_reported_1y = zeros(numel(age_vec), 2);
num_para_sets = size(ascertainment_ratio, 1);
num_ages = numel(age_vec);
num_stages = mc.num_CRC_stages;

ascertainment_ratio_this = zeros(num_para_sets, num_ages, num_stages);

age_knots = age_range(1) + diff(age_range) * ascertainment_ratio(:, 1:3:end);

for CRC_stage = 1:num_stages

    idx = (CRC_stage - 1) * 3 + (1:3);
    asr = ascertainment_ratio(:, idx);

    asr_low = asr(:, 2);
    asr_high = asr(:, 3);

    is_low = age_vec < age_knots(:, CRC_stage);
    is_high = age_vec > age_range(end);
    is_mid = ~is_low & ~is_high;

    ascertainment_ratio_this(:, is_low, CRC_stage) = repmat(asr_low, 1, sum(is_low));

    ascertainment_ratio_this(:, is_high, CRC_stage) = repmat(asr_high, 1, sum(is_high));

    mid_ages = age_vec(is_mid);
    mid_vals = pchip([age_knots(:, CRC_stage), age_range(end)], [asr_low, asr_high]', mid_ages);
    ascertainment_ratio_this(:, is_mid, CRC_stage) = mid_vals';
end

cohort.ascertainment_ratio_this=ascertainment_ratio_this;
cohort.incidence_reported_1y = zeros(numel(age_vec), num_CRC_stages+1);

for CRC_stage = 1:num_CRC_stages
    incidence_this = cohort.incidence_CRC_1y(:, CRC_stage);
    cohort.incidence_reported_1y(:, CRC_stage) = ascertainment_ratio_this(:,:,CRC_stage)'.*incidence_this;
    cohort.incidence_reported_1y(:, num_CRC_stages+1) = cohort.incidence_reported_1y(:, num_CRC_stages+1)+(1-ascertainment_ratio_this(:,:, CRC_stage)').*incidence_this;
end

cohort.age_5y = mc.age_knots_HKCR;
for ii = 2:numel(cohort.age_5y)
    fI = find(cohort.age_1y >= cohort.age_5y(ii-1) & cohort.age_1y <= cohort.age_5y(ii));
    if (isempty(fI))
        continue;
    end
    cohort.size_5y(ii-1, :) = sum(cohort.size_1y(fI, :), 1)-sum(cohort.size_1y(fI([1 end]), :), 1)/2;
    cohort.incidence_reported_5y(ii-1, :) = sum(cohort.incidence_reported_1y((ii-1)*5+1:min(end, 5*ii), :));
    cohort.stats_5y{CRC_INCIDENCE, MODEL}(ii-1, :) = cohort.incidence_reported_5y(ii-1, :);
    cohort.stats_5y{CRC_INCIDENCE_RATE, MODEL}(ii-1, :) = cohort.incidence_reported_5y(ii-1, :)./repmat(cohort.size_5y(ii-1, COHORT_TOTAL), 1, num_CRC_stages+1);

    CRC_mortality_5y = cohort.state_1y(fI(end), CRC_DEATH)-cohort.state_1y(fI(1), CRC_DEATH);
    cohort.stats_5y{CRC_MORTALITY, MODEL}(ii-1) = CRC_mortality_5y;
    cohort.stats_5y{CRC_MORTALITY_RATE, MODEL}(ii-1) = CRC_mortality_5y/cohort.size_5y(ii-1, COHORT_TOTAL);
end

cohort.prevalence_UCRC_1y=cohort.prevalence_1y(:,:,COHORT_WITHOUT_DIAGNOSED_CRC);
cohort.prevalence_DC_1y=cohort.prevalence_1y(:,:,COHORT_TOTAL);

end
