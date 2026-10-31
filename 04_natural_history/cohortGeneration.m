function dydt = cohortGeneration(tt, yy,age_knots_HKCR,all_cause_mortality,ini,mor)

all_cause_mortality_rate = interp1(age_knots_HKCR,all_cause_mortality, tt);

for jj=1:size(ini, 2)

    CRC_incidence_rate(jj)= interp1(age_knots_HKCR,ini(:,jj), tt);

end

CRC_mortality_rate = interp1(age_knots_HKCR, mor, tt);

dydt = [
    -all_cause_mortality_rate
    CRC_incidence_rate'
    CRC_mortality_rate
    all_cause_mortality_rate
]*yy(1);

end