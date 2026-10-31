function inm=prepareObservedCohort(mc,gender_this)

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

all_cause_mortality=cohort(gender_this).interpolated_rate{mc.ALL_CAUSE_MORTALITY};
CRC_incidence_rate=cohort(gender_this).interpolated_rate{mc.CRC_INCIDENCE};
CRC_mortality_rate=cohort(gender_this).interpolated_rate{mc.CRC_MORTALITY};

[tSln, xSln] = ode45(@(tt,yy)cohortGeneration(tt, yy,mc.age_knots_HKCR,all_cause_mortality,CRC_incidence_rate,CRC_mortality_rate), ...
   mc.age_range, [mc.cohort_size 0 0 0 0 0 0 0 0]);

    one_year_age_groups = mc.age_range(1):mc.age_range(end);
    ten_year_age_groups = mc.age_knots_HKCR;
    ys = interp1(tSln, xSln, one_year_age_groups);
    inm.obs_cohort_size = ys;

    for ii = 2:numel(ten_year_age_groups)
        fI = find(one_year_age_groups >= ten_year_age_groups(ii-1) & one_year_age_groups <= ten_year_age_groups(ii));
        group_size = sum(ys(fI, 1))-sum(ys(fI([1 end])))/2;
        for jj=1:size(CRC_incidence_rate, 2)
        incidence_this = ys(fI(end), 1+jj)-ys(fI(1), 1+jj);
        inm.stats_5y{mc.CRC_INCIDENCE, mc.OBSERVED}(ii-1,jj) = incidence_this;
        inm.stats_5y{mc.CRC_INCIDENCE_RATE, mc.OBSERVED}(ii-1,jj) = incidence_this/group_size;
        end
        CRC_mortality_this = ys(fI(end), end-1)-ys(fI(1), end-1);
        inm.stats_5y{mc.CRC_MORTALITY, mc.OBSERVED}(ii-1) = CRC_mortality_this;
        inm.stats_5y{mc.CRC_MORTALITY_RATE, mc.OBSERVED}(ii-1) = CRC_mortality_this/group_size;
    end

end