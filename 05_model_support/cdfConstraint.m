function fun_val = cdfConstraint(model_specs, mc)

fun_val = 0;
cdf_constraint = mc.cdf_constraint;

NORMAL = mc.NORMAL;
LSL = mc.LSL;
NA = mc.NA;
HSL = mc.HSL;
AA = mc.AA;
CRC = mc.CRC;

vec_v = linspace(1e-100, 10, 100);

CRC_pathway{1} = [NORMAL NA AA CRC];
CRC_pathway{2} = [NORMAL LSL HSL CRC];

check_ages = [30, 50, 70];

if (cdf_constraint >= 1)
    for age_check = check_ages
        for contractor_idx = mc.contractor_list
            for pathway_idx = 1:numel(CRC_pathway)
                pathway_this = CRC_pathway{pathway_idx};
                for stage_idx = 2:numel(pathway_this)
                    stage_2 = pathway_this(stage_idx);
                    stage_1 = pathway_this(stage_idx-1);
                    if (stage_1==NORMAL)
                        cdf_obs{contractor_idx}(stage_1, :) = cdfUnderlyingLogFIT(vec_v, stage_1, contractor_idx, model_specs, age_check);
                    end
                    cdf_obs{contractor_idx}(stage_2, :) = cdfUnderlyingLogFIT(vec_v, stage_2, contractor_idx, model_specs, age_check);
                    delta = 1-cdf_obs{contractor_idx}(stage_1, :)./cdf_obs{contractor_idx}(stage_2, :);
                    if (any(delta>1e-100, 'all'))
                        fun_val = -100000;
                        return;
                    end
                end
            end
        end
    end
end

if (cdf_constraint >= 2)
    for age_check = check_ages
        for contractor_idx = mc.contractor_list
            for pathway_idx = 1:numel(CRC_pathway)
                pathway_this = CRC_pathway{pathway_idx};
                for stage_idx = 2:numel(pathway_this)
                    stage_2 = pathway_this(stage_idx);
                    stage_1 = pathway_this(stage_idx-1);
                    if (stage_1==NORMAL)
                        cdf_obs{contractor_idx}(stage_1, :) = cdfObservedLogFIT(vec_v, stage_1, contractor_idx, model_specs, age_check);
                    end
                    cdf_obs{contractor_idx}(stage_2, :) = cdfObservedLogFIT(vec_v, stage_2, contractor_idx, model_specs, age_check);
                    delta = 1-cdf_obs{contractor_idx}(stage_1, :)./cdf_obs{contractor_idx}(stage_2, :);
                    if (any(delta>1e-100, 'all'))
                        fun_val = -100000;
                        return;
                    end
                end
            end
        end
    end
end

if (cdf_constraint >= 3)
    for age_check = check_ages
        for contractor_idx = mc.contractor_list
            for pathway_idx = 1:numel(CRC_pathway)
                pathway_this = CRC_pathway{pathway_idx};
                for stage_idx = 2:numel(pathway_this)
                    stage_2 = pathway_this(stage_idx);
                    stage_1 = pathway_this(stage_idx-1);
                    if (stage_1==NORMAL)
                        cdf_obs_x2_max{contractor_idx}(stage_1, :) = cdfMaxTwoObservedLogFIT(vec_v, stage_1, contractor_idx, model_specs, age_check);
                   end
                    cdf_obs_x2_max{contractor_idx}(stage_2, :) = cdfMaxTwoObservedLogFIT(vec_v, stage_2, contractor_idx, model_specs, age_check);
                    delta = 1-cdf_obs_x2_max{contractor_idx}(stage_1, :)./cdf_obs_x2_max{contractor_idx}(stage_2, :);
                    if (any(delta>1e-100, 'all'))
                        fun_val = -100000;
                        return;
                    end
                end
            end
        end
    end
end
