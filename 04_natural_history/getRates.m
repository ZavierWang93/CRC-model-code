function cohort=getRates(gender_this,mc)
    for gender = gender_this

        if gender==1
        cohort(gender).interpolated_rate{mc.ALL_CAUSE_MORTALITY}(1) = 0.3/1e3;
        else
        cohort(gender).interpolated_rate{mc.ALL_CAUSE_MORTALITY}(1) = 0.15/1e3;
        end

        for jj = 1:size(mc.CRC_incidence_rate_HKCR{1}, 2)
            cohort(gender).interpolated_rate{mc.CRC_INCIDENCE}(1, jj) = 0;
        end

        cohort(gender).interpolated_rate{mc.CRC_MORTALITY}(1) = 0;
        options = optimset('TolX', 1e-16);

        for ii = 2:numel(mc.age_knots_HKCR)
            para.t0 = mc.age_knots_HKCR(ii-1);
            para.t1 = mc.age_knots_HKCR(ii);

            para.a0 = cohort(gender).interpolated_rate{mc.ALL_CAUSE_MORTALITY}(ii-1);
            para.mRate = mc.all_cause_mortality_rate_HKCSD(ii-1, gender);
            para.a1 = fminbnd(@(x)cohortObjective1(x, para), para.a0, 50*para.mRate, options);
            cohort(gender).interpolated_rate{mc.ALL_CAUSE_MORTALITY}(ii) = para.a1;

            for jj = 1:size(mc.CRC_incidence_rate_HKCR{1}, 2)
                para.b0 = cohort(gender).interpolated_rate{mc.CRC_INCIDENCE}(ii-1, jj);
                para.incRate = mc.CRC_incidence_rate_HKCR{gender}(ii-1, jj);
                para.b1 = fminbnd(@(x)cohortObjective2(x, para), 1e-10, 50*para.incRate, options);
                cohort(gender).interpolated_rate{mc.CRC_INCIDENCE}(ii, jj) = para.b1;
            end

            para.b0 = cohort(gender).interpolated_rate{mc.CRC_MORTALITY}(ii-1);
            para.incRate = mc.CRC_mortality_rate_HKCR(ii-1, gender);
            para.b1 = fminbnd(@(x)cohortObjective2(x, para), 1e-10, 50*para.incRate, options);
            cohort(gender).interpolated_rate{mc.CRC_MORTALITY}(ii) = para.b1;
        end
    end
end