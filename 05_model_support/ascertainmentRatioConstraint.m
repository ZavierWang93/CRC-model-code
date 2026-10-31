function fun_val = ascertainmentRatioConstraint(model_specs, mc)

    age_range = [20, 100];
    age_vec = age_range(1):age_range(end);

    ascertainment_ratio = [model_specs.ascertainment_ratio model_specs.ascertainment_ratio(7:9)];
    num_para_sets = size(ascertainment_ratio, 1);
    num_ages = length(age_vec);
    num_stages = mc.num_CRC_stages;

    ascertainment_ratio_this = zeros(num_para_sets, num_ages, num_stages);

    age_knots = age_range(1) + diff(age_range) * ascertainment_ratio(:, 1:3:end);

    for CRC_stage = 1:num_stages

        idx = (CRC_stage - 1) * 3 + (1:3);
        asr = ascertainment_ratio(:, idx);

        age_knot = age_knots(:, CRC_stage);
        asr_low = asr(:, 2);
        asr_high = asr(:, 3);

        is_low = age_vec < age_knot';
        is_high = age_vec > age_range(end);
        is_mid = ~is_low & ~is_high;

        ascertainment_ratio_this(:, is_low, CRC_stage) = repmat(asr_low, 1, sum(is_low));

        ascertainment_ratio_this(:, is_high, CRC_stage) = repmat(asr_high, 1, sum(is_high));

        mid_ages = age_vec(is_mid);
        mid_vals = pchip([age_knot, age_range(end)], [asr_low, asr_high]', mid_ages);
        ascertainment_ratio_this(:, is_mid, CRC_stage) = mid_vals';
    end

    delta = diff(ascertainment_ratio_this, 1, 2);
    if any(any(delta > 0))
        fun_val1 = -1e5;
    else
        fun_val1 = 0;
    end

    isValid = (ascertainment_ratio_this(:,:,4) >= ascertainment_ratio_this(:,:,3)) &...
        (ascertainment_ratio_this(:,:,3) >= ascertainment_ratio_this(:,:,2)) &...
        (ascertainment_ratio_this(:,:,2) >= ascertainment_ratio_this(:,:,1));

    if any(any(isValid == 0 ))
        fun_val2 = -1e5;
    else
        fun_val2 = 0;
    end

    fun_val=fun_val1+fun_val2;

end
