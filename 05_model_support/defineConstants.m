NORMAL = 1;
LSL = 2;
NA = 3;
HSL = 4;
AA = 5;
CRC = 6;
UNKNOWN = 7;

main_control.NORMAL = NORMAL;
main_control.LSL = LSL;
main_control.NA = NA;
main_control.HSL = HSL;
main_control.AA = AA;
main_control.CRC = CRC;
main_control.UNKNOWN = UNKNOWN;
main_control.num_outcomes = UNKNOWN;
main_control.num_stages = CRC;

main_control.outcome_tag = {'N', 'LSL', 'NA', 'HSL', 'AA', 'CRC', 'U'};

MALE = 1;
FEMALE = 2;
main_control.MALE = MALE;
main_control.FEMALE = FEMALE;
main_control.gender_tag = {'Male', 'Female'};

num_CRC_stages = 4;
num_CRC_diagnosed_substages = 1;
death_type = 2;

main_control.num_CRC_stages = num_CRC_stages;
main_control.num_CRC_diagnosed_substages = num_CRC_diagnosed_substages;
main_control.death_type = death_type;

CRC_UNDIAGNOSED = AA + (1:num_CRC_stages);
if (main_control.death_type == 1)
    CRC_DIAGNOSED = CRC_UNDIAGNOSED(end) + (1:num_CRC_diagnosed_substages * num_CRC_stages);
else
    CRC_DIAGNOSED = CRC_UNDIAGNOSED(end) + (1:num_CRC_stages);
end
CRC_SURVIVOR = CRC_DIAGNOSED(end) + 1;
CRC_DEATH = CRC_SURVIVOR + 1;
DEATH = CRC_DEATH + 1;
CUM_CRC_DIAGNOSED = DEATH + (1:num_CRC_stages);
CRC_AA = CUM_CRC_DIAGNOSED(end) + 1;
CRC_SL = CRC_AA + 1;

main_control.CRC_UNDIAGNOSED = CRC_UNDIAGNOSED;
main_control.CRC_DIAGNOSED = CRC_DIAGNOSED;
main_control.CRC_SURVIVOR = CRC_SURVIVOR;
main_control.CRC_DEATH = CRC_DEATH;
main_control.DEATH = DEATH;
main_control.CUM_CRC_DIAGNOSED = CUM_CRC_DIAGNOSED;
main_control.CRC_AA = CRC_AA;
main_control.CRC_SL = CRC_SL;

main_control.num_exit_rates = numel([NA AA LSL HSL CRC_UNDIAGNOSED]);

main_control.cohort_size = 100000;
main_control.age_range = [20 100];
main_control.age_knots_HKCR = [20 30:5:85 100];
main_control.MODEL = 2;
main_control.OBSERVED = 1;

main_control.COHORT_TOTAL = 1;
main_control.COHORT_WITHOUT_DIAGNOSED_CRC = 2;

main_control.CRC_INCIDENCE = 1;
main_control.CRC_MORTALITY = 2;
main_control.CRC_INCIDENCE_RATE = 3;
main_control.CRC_MORTALITY_RATE = 4;
main_control.ALL_CAUSE_MORTALITY = 5;
main_control.ALL_CAUSE_MORTALITY_RATE = 6;

main_control.contractor_list = [1 2];
main_control.knot_incidence_rate = [20 45 70];
main_control.knot_age_effect_disease_progression = [50 75 100];
main_control.knot_ascertainment_ratio = [50 100];

INVALID_FIT = -0.5;
main_control.INVALID_FIT = INVALID_FIT;
main_control.ZERO = 1;
main_control.EXACT = 2;
main_control.RIGHT_CENSORED = 3;
main_control.LEFT_CENSORED = 4;
main_control.INVALID = 5;

FIT_right_censored = [120000 1000];
FIT_left_censored = [INVALID_FIT/2 9];
FIT_pos_threshold = 100;

main_control.logFIT_pos_threshold = log(1 + FIT_pos_threshold);
main_control.logFIT_right_censored = log(1 + FIT_right_censored);
main_control.logFIT_left_censored = log(1 + FIT_left_censored);

CRC_incidence_rate_HKCR_raw = [
    0.1 	0.7 	0.3 	2.1 	2.9 	6.1 	12.0 	15.2 	27.9 	27.9 	43.4 	37.8 	22.3
    0.1 	0.4 	0.4 	1.6 	2.9 	5.7 	8.0 	12.5 	18.6 	23.9 	22.0 	25.6 	15.2
    0.1 	0.5 	0.4 	1.8 	2.9 	5.9 	10.0 	13.8 	23.4 	25.9 	32.2 	30.9 	17.5
    0.2 	0.9 	2.7 	3.5 	7.3 	14.6 	22.9 	40.4 	55.4 	79.2 	111.3 	123.1 	96.3
    0.1 	0.6 	1.1 	2.9 	4.0 	8.5 	17.7 	24.1 	32.0 	50.1 	74.4 	82.5 	50.3
    0.1 	0.7 	1.8 	3.2 	5.5 	11.5 	20.3 	32.2 	43.9 	64.9 	91.9 	100.1 	65.3
    0.3 	1.0 	1.6 	5.7 	10.2 	18.1 	30.0 	42.1 	67.6 	78.0 	102.6 	106.6 	66.3
    0.2 	0.8 	2.2 	4.9 	7.8 	14.4 	21.7 	28.8 	36.7 	52.9 	64.5 	80.6 	46.4
    0.2 	0.9 	2.0 	5.2 	8.8 	16.2 	25.8 	35.4 	52.5 	65.7 	82.6 	91.9 	52.9
    0.3 	1.2 	2.0 	4.0 	8.1 	14.3 	27.7 	40.8 	56.9 	73.9 	92.1 	113.3 	114.0
    0.1 	0.6 	2.2 	4.1 	6.6 	11.8 	17.8 	21.2 	33.2 	37.5 	62.4 	68.3 	65.9
    0.2 	0.9 	2.1 	4.1 	7.3 	13.0 	22.7 	30.9 	45.3 	56.0 	76.5 	87.8 	81.6
    0.3 	2.1 	2.8 	5.4 	6.9 	12.0 	14.9 	21.7 	29.8 	43.8 	72.2 	111.2 	163.0
    0.2 	0.7 	1.8 	3.6 	5.2 	6.3 	9.7 	11.3 	19.0 	21.9 	38.6 	69.2 	125.5
    0.2 	1.3 	2.2 	4.3 	6.0 	9.0 	12.2 	16.5 	24.5 	33.1 	54.5 	87.4 	137.8
    1.1 	5.9 	9.5 	20.7 	35.3 	65.1 	107.5 	160.2 	237.6 	302.9 	421.6 	491.9 	461.9
    0.6 	3.1 	7.8 	17.0 	26.5 	46.7 	74.9 	97.9 	139.5 	186.2 	261.8 	326.3 	303.4
    0.9 	4.2 	8.5 	18.6 	30.5 	55.6 	91.1 	128.9 	189.6 	245.6 	337.8 	398.1 	355.1
]/1e5;

main_control.CRC_incidence_rate_HKCR{MALE} = CRC_incidence_rate_HKCR_raw(1:3:end, :)';
main_control.CRC_incidence_rate_HKCR{FEMALE} = CRC_incidence_rate_HKCR_raw(2:3:end, :)';

main_control.CRC_mortality_rate_HKCR = [
    0.6 0.1
    1.6	0.9
    2.3	1.5
    4.2	3.8
    10.1	6.8
    16.6	12.8
    32.5	20.8
    53.2	29.7
    70.7	43.5
    126.1	64.3
    190.4	111.1
    269.6	175.2
    407.8	272
]/1e5;

main_control.all_cause_mortality_rate_HKCSD = [
    0.4	0.2
    0.6	0.3
    0.9	0.45
    1.3	0.7
    2.05	1.05
    3.3	1.8
    5.35	2.6
    8.5	4.05
    13.85	5.9
    22	10.15
    36.8	18.8
    64.7	35.8
    126.3	92.6
]/1e3;

main_control.CRC_5y_survival_prob = [0.96 0.87 0.69 0.093];

colonoscopy_sensitivity = zeros(1, CRC);
colonoscopy_sensitivity(NORMAL) = 1;
colonoscopy_sensitivity(LSL) = 0.7;
colonoscopy_sensitivity(NA) = 0.8;
colonoscopy_sensitivity(HSL) = 0.8;
colonoscopy_sensitivity(AA) = 0.9;
colonoscopy_sensitivity(CRC) = 0.95;
main_control.colonoscopy_sensitivity = colonoscopy_sensitivity;
