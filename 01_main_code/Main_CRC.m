clearvars;
close all;
rng('shuffle');

project_root = 'D:\Dropbox\CRC SCREENING PROJECT\Manuscript Submit\Base_Model_Code_20261029';
data_root = 'D:\Dropbox\CRC SCREENING PROJECT\Manuscript Submit\Base_Model_Code_20261029';
addpath(genpath(project_root));
chain_id = 1;
resume_run = false;
Step_size = [];

if isempty(gcp('nocreate')) || ~isa(gcp('nocreate'), 'parallel.ThreadPool')
    delete(gcp('nocreate'));
    parpool('Threads');
end

defineConstants;

main_control.preclinical_constraint = 0;
main_control.cdf_constraint= 2;
main_control.SL_constraint = 2;
main_control.Prevalence_constraint = 1;
main_control.No_CS_pathway_pro=0.95;

for stage=NORMAL:CRC
    main_control.logFIT_dist{stage} = 'Lognormal';
end

prevalence = [1 1e-3 1e-3 1e-3 1e-3 1e-3];
prevalence = prevalence/sum(prevalence);
main_control.num_prevalence=numel(prevalence);

ascertainment_ratio=repmat(0.8,1,9); 
main_control.num_ascertainment_ratio=numel(ascertainment_ratio);

inc_NA = pchip([20 75], [1e-4, 0.1], main_control.knot_incidence_rate);
main_control.num_inc_NA=numel(inc_NA);

inc_LSL = pchip([20 75], [1e-4, 0.03], main_control.knot_incidence_rate);
main_control.num_inc_LSL=numel(inc_LSL);

incidence_rate = [inc_NA inc_LSL];

exit_rate = repmat(0.1, 1, numel([LSL NA HSL AA CRC_UNDIAGNOSED])); 
exit_rate(end-main_control.num_CRC_stages+1:end) =  [0.5 0.5 0.5 0.5];
main_control.num_exit_rate=numel(exit_rate);

age_effect_on_disease_progression = repmat(3, 1, numel(main_control.knot_age_effect_disease_progression)-1)./10;
main_control.num_age_effect_on_disease_progression=numel(age_effect_on_disease_progression);

prob_diagnosed = [0.1 0.3 0.5];
main_control.num_prob_diagnosed=numel(prob_diagnosed);

knot_using1=(linspace(log(10), log(1000), 5));
knot_omega{1} =exp([log(10),knot_using1(2:4),log(1000)]);

knot_using2=(linspace(log(10), log(1000), 5));
knot_omega{2} = exp([log(10),knot_using2(2:4),log(1000)]);

omega{1} = [0.5 1 2 1 1]./10;
omega{2} = [0.5 1 2 1 1]./10;
main_control.num_omega{1}=numel(omega{1});
main_control.num_omega{2}=numel(omega{2});

main_control.num_knot_omega{1} = numel(knot_omega{1});
main_control.num_knot_omega{2} = numel(knot_omega{2});

mean_v = [3 3.5 4 4.5 5 8]./10;
main_control.num_mean_v=numel(mean_v);
std_v = [0.5 0.8 1.1 1.4 1.7 2]./10;
main_control.num_std_v=numel(std_v);
pZero = [0.1 0.08 0.06 0.04 0.02 0.01];
main_control.num_pZero=numel(pZero);
q0 = 0.01;
main_control.num_q0=numel(q0);
qNA = 0.03;
main_control.num_qNA=numel(qNA);
age_eff_mean = repmat(1e-6, 1, numel(mean_v))./10;
main_control.num_age_eff_mean = numel(age_eff_mean);
main_control.knot_omega = knot_omega;

load('D:\Dropbox\CRC SCREENING PROJECT\dataset.mat');

for choosing_gender=1:2
inm{choosing_gender} =prepareObservedCohort(main_control,gender_this);
end

pathway = fullfile(data_root, 'Final Results', ...
    sprintf('mcmc_base_Chain%d.mat', chain_id));
main_control.continuous_scenario_id = 'base';
main_control.continuous_data_id = 'original_path_data_and_filters';
main_control.resume = resume_run;
main_control.max_production_samples = Inf;
main_control.skip_tuning = false;
main_control.tuning_min_cycles = 1000;
main_control.tuning_max_cycles = 10000;
main_control.tuning_window = 200;
main_control.save_interval_sec = 600;
main_control.chunk_size = 1000;
main_control.tail_size = 10000;
main_control.report_interval = 200;
main_control.diagnostic_interval = 500;
main_control.diagnostic_min_samples = 1000;
main_control.ESS_target = 400;
comparison_chain_id = 1;
if chain_id == 1, comparison_chain_id = 2; end
other_pathway = strrep(pathway, sprintf('_Chain%d.mat',chain_id), ...
    sprintf('_Chain%d.mat',comparison_chain_id));
main_control.convergence_paths = {pathway,other_pathway};

%%
tic
answer = MCMC_CRC( ...
    Step_size, main_control, ini_1, ini_2, ...
    sorted_obs_logFIT{1}, sorted_obs_logFIT{2}, inm{1}, inm{2}, ...
    pathway, likelihood_fun);
toc
