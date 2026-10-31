function answer = MCMCContinuous_CRC( ...
        Step_size, mc, initial_para_M, initial_para_F, ...
        observe_sorts_M, observe_sorts_F, inm_M, inm_F, pathway, lh_fun)

if nargin < 10 || isempty(lh_fun)
    error('MCMC_meanage_continuous:MissingLikelihood', ...
        'An explicit meanage likelihood function is required.');
end
if ischar(lh_fun) || isstring(lh_fun)
    lh_fun = str2func(lh_fun);
end
if ~isa(lh_fun, 'function_handle')
    error('MCMC_meanage_continuous:InvalidLikelihood', ...
        'lh_fun must be a function handle or function name.');
end
if nargin < 2 || isempty(mc)
    mc = struct();
end

initial_para_M = double(initial_para_M(:)');
initial_para_F = double(initial_para_F(:)');
nP = numel(initial_para_M);
if numel(initial_para_F) ~= nP || ~ismember(nP, [64, 99, 105, 106, 113])
    error('MCMC_meanage_continuous:InvalidParameterCount', ...
        'Male and female vectors must have the same supported length: 64, 99, 105, 106 or 113.');
end

layout = mcmcBlockLayout(nP, mc);
settings = readSettings(mc, Step_size, layout);
paths = buildPaths(pathway);
lockCleanup = acquireLock(paths.lock, settings.force_unlock);
identity = buildIdentity(mc, layout, lh_fun, settings);

if settings.resume
    runtime = loadCheckpoint(paths.checkpoint, identity, settings);
    runtime.state.mc.convergence_paths = controlValue(mc, 'convergence_paths', {});
    rng(runtime.rng_state);
    fprintf('Resuming %s at production sample %d (%s).\n', ...
        identity.scenario_id, runtime.total_production_samples, runtime.stage);
else
    if isfile(paths.checkpoint) || isfile(paths.snapshot)
        error('MCMC_meanage_continuous:OutputExists', ...
            'Output already exists. Use a new pathway or set mc.resume=true: %s', ...
            paths.snapshot);
    end
    runtime = initializeRuntime(initial_para_M, initial_para_F, mc, ...
        observe_sorts_M, observe_sorts_F, inm_M, inm_F, layout, settings, ...
        identity, lh_fun);
    commitState(runtime, paths, settings);
end

if strcmp(runtime.stage, 'tuning')
    runtime = runTuning(runtime, observe_sorts_M, observe_sorts_F, ...
        inm_M, inm_F, layout, settings, identity, lh_fun, paths);
    if strcmp(runtime.status, 'tuning_paused')
        answer = buildAnswer(runtime, paths, settings);
        return
    end
end

if strcmp(runtime.status, 'paused')
    runtime.status = 'production_running';
    runtime.stage = 'production';
end
runtime.stage = 'production';
runtime.status = 'production_running';
lastCommit = tic;

fprintf(['Continuous meanage MCMC: %d parameters/sex, %d active blocks, ' ...
    '%d proposals per stored cycle.\n'], nP, nnz(layout.active_mask), ...
    layout.proposals_per_macro);
fprintf('  Proposal: independent uniform random walk on the original scale.\n');
fprintf('  Likelihood: %s\n', identity.likelihood_name);

while isinf(settings.max_production_samples) || ...
        runtime.total_production_samples < settings.max_production_samples
    cycleTimer = tic;
    [runtime.state, accepted, attempted, decisions] = runMacroCycle( ...
        runtime.state, layout, runtime.step_size, lh_fun, ...
        observe_sorts_M, observe_sorts_F, inm_M, inm_F);
    elapsed = toc(cycleTimer);
    runtime = appendProduction(runtime, accepted, attempted, decisions, elapsed, ...
        settings, nP);

    if mod(runtime.total_production_samples, settings.report_interval) == 0
        runtime = reportAcceptance(runtime, layout);
    end
    if runtime.total_production_samples >= settings.diagnostic_min_samples && ...
            mod(runtime.total_production_samples, settings.diagnostic_interval) == 0
        runtime = updateDiagnostics(runtime, layout, settings);
    end

    if runtime.pending.count == settings.chunk_size
        runtime = flushPending(runtime, paths, settings);
        commitState(runtime, paths, settings);
        lastCommit = tic;
    elseif toc(lastCommit) >= settings.save_interval_sec
        runtime.rng_state = rng;
        commitState(runtime, paths, settings);
        lastCommit = tic;
    end

    if isfile(paths.pause)
        runtime = flushPending(runtime, paths, settings);
        runtime.status = 'paused';
        runtime.stage = 'production';
        runtime.rng_state = rng;
        commitState(runtime, paths, settings);
        delete(paths.pause);
        fprintf('Paused safely after production sample %d.\n', ...
            runtime.total_production_samples);
        answer = buildAnswer(runtime, paths, settings);
        return
    end
end

runtime = flushPending(runtime, paths, settings);
runtime.status = 'production_complete';
runtime.stage = 'production';
runtime.rng_state = rng;
commitState(runtime, paths, settings);
answer = buildAnswer(runtime, paths, settings);
end

function settings = readSettings(mc, Step_size, layout)
settings = struct();
settings.max_production_samples = controlValue(mc, ...
    'max_production_samples', Inf);
numBlocks = numel(layout.blocks);
settings.skip_tuning = logicalControl(mc, 'skip_tuning', false);
settings.resume = logicalControl(mc, 'resume', false);
settings.force_unlock = logicalControl(mc, 'force_unlock', false);
settings.tuning_window = controlValue(mc, 'tuning_window', 200);
settings.tuning_min_cycles = controlValue(mc, 'tuning_min_cycles', 400);
settings.tuning_max_cycles = controlValue(mc, 'tuning_max_cycles', 2000);
settings.accept_min = controlValue(mc, 'accept_min', 0.20);
settings.accept_max = controlValue(mc, 'accept_max', 0.70);
settings.scale_down = controlValue(mc, 'scale_down', 0.60);
settings.scale_up = controlValue(mc, 'scale_up', 1.40);
settings.chunk_size = controlValue(mc, 'chunk_size', 1000);
settings.tail_size = controlValue(mc, 'tail_size', 10000);
settings.save_interval_sec = controlValue(mc, 'save_interval_sec', 900);
settings.report_interval = controlValue(mc, 'report_interval', 200);
settings.diagnostic_interval = controlValue(mc, 'diagnostic_interval', 500);
settings.diagnostic_min_samples = controlValue(mc, ...
    'diagnostic_min_samples', 1000);
settings.ess_target = controlValue(mc, 'ESS_target', 400);
settings.ess_stop_pct = controlValue(mc, 'ESS_stop_pct', 0.95);
settings.max_diagnostic_records = controlValue(mc, ...
    'max_diagnostic_records', 200);
settings.tuning_block_mask = controlValue(mc, ...
    'tuning_block_mask', layout.active_mask);

integerFields = {'tuning_window', 'tuning_min_cycles', ...
    'tuning_max_cycles', 'chunk_size', 'tail_size', 'report_interval', ...
    'diagnostic_interval', 'diagnostic_min_samples', ...
    'max_diagnostic_records'};
for k = 1:numel(integerFields)
    value = settings.(integerFields{k});
    if ~isscalar(value) || ~isfinite(value) || value < 1 || value ~= round(value)
        error('MCMC_meanage_continuous:InvalidSetting', ...
            'mc.%s must be a positive integer.', integerFields{k});
    end
end
if ~(isscalar(settings.max_production_samples) && (settings.max_production_samples == Inf || ...
        (isscalar(settings.max_production_samples) && ...
         isfinite(settings.max_production_samples) && ...
         settings.max_production_samples >= 0 && ...
         settings.max_production_samples == round(settings.max_production_samples))))
    error('MCMC_meanage_continuous:InvalidProductionLimit', ...
        'mc.max_production_samples must be Inf or a nonnegative integer.');
end
if settings.tuning_min_cycles > settings.tuning_max_cycles || ...
        mod(settings.tuning_min_cycles, settings.tuning_window) ~= 0 || ...
        mod(settings.tuning_max_cycles, settings.tuning_window) ~= 0
    error('MCMC_meanage_continuous:InvalidTuningSchedule', ...
        ['Tuning limits must be multiples of tuning_window, and the ' ...
         'minimum cannot exceed the maximum.']);
end
if settings.accept_min < 0 || settings.accept_max > 1 || ...
        settings.accept_min >= settings.accept_max
    error('MCMC_meanage_continuous:InvalidAcceptanceRange', ...
        'Acceptance limits must satisfy 0 <= accept_min < accept_max <= 1.');
end
if ~isscalar(settings.scale_down) || ~isfinite(settings.scale_down) || ...
        settings.scale_down <= 0 || settings.scale_down >= 1 || ...
        ~isscalar(settings.scale_up) || ~isfinite(settings.scale_up) || ...
        settings.scale_up <= 1
    error('MCMC_meanage_continuous:InvalidScaleFactor', ...
        'scale_down must be in (0,1), and scale_up must exceed 1.');
end
if ~isscalar(settings.save_interval_sec) || ...
        ~(isinf(settings.save_interval_sec) || ...
          (isfinite(settings.save_interval_sec) && settings.save_interval_sec > 0))
    error('MCMC_meanage_continuous:InvalidSaveInterval', ...
        'save_interval_sec must be positive or Inf.');
end
if ~isscalar(settings.ess_target) || ~isfinite(settings.ess_target) || ...
        settings.ess_target <= 0 || ...
        ~isscalar(settings.ess_stop_pct) || ~isfinite(settings.ess_stop_pct) || ...
        settings.ess_stop_pct < 0 || settings.ess_stop_pct > 1
    error('MCMC_meanage_continuous:InvalidESSSetting', ...
        'ESS_target must be positive and ESS_stop_pct must lie in [0,1].');
end
if settings.tail_size < settings.chunk_size
    error('MCMC_meanage_continuous:InvalidTailSize', ...
        'tail_size must be at least chunk_size.');
end
if ~islogical(settings.tuning_block_mask)
    if ~isnumeric(settings.tuning_block_mask) || ...
            any(~ismember(settings.tuning_block_mask(:), [0, 1]))
        error('MCMC_meanage_continuous:InvalidTuningMask', ...
            'mc.tuning_block_mask must contain one logical or 0/1 value per block.');
    end
    settings.tuning_block_mask = logical(settings.tuning_block_mask);
end
if ~isvector(settings.tuning_block_mask) || ...
        numel(settings.tuning_block_mask) ~= numBlocks
    error('MCMC_meanage_continuous:InvalidTuningMask', ...
        'mc.tuning_block_mask must contain one value per canonical block.');
end
settings.tuning_block_mask = reshape(settings.tuning_block_mask, 1, []);
if any(settings.tuning_block_mask & ~layout.active_mask)
    error('MCMC_meanage_continuous:InvalidTuningMask', ...
        'Inactive blocks cannot be selected for tuning.');
end
if ~any(settings.tuning_block_mask)
    settings.skip_tuning = true;
end

if isempty(Step_size)
    settings.initial_step_size = 0.01*ones(1, numBlocks);
    settings.initial_step_size(strcmp(layout.block_names, 'age_eff_mean')) = 0.0002;
elseif isscalar(Step_size)
    settings.initial_step_size = repmat(double(Step_size), 1, numBlocks);
elseif isvector(Step_size) && numel(Step_size) == numBlocks
    settings.initial_step_size = double(Step_size(:)');
elseif layout.parameter_count == 105 && numel(Step_size) == 19

    legacy = double(Step_size(:)');
    settings.initial_step_size = [legacy(1:4) legacy(4) legacy(5:end)];
else
    error('MCMC_meanage_continuous:InvalidStepSize', ...
        'Step_size must be empty, scalar, or match answer.block_names.');
end
activeSteps = settings.initial_step_size(layout.active_mask);
if any(~isfinite(activeSteps) | activeSteps <= 0)
    error('MCMC_meanage_continuous:InvalidStepSize', ...
        'Every active block step size must be finite and positive.');
end
end

function paths = buildPaths(pathway)
if isempty(pathway)
    error('MCMC_meanage_continuous:MissingPathway', ...
        'A snapshot pathway is required for continuous sampling.');
end
pathway = char(string(pathway));
[parent, name, ext] = fileparts(pathway);
if isempty(ext)
    pathway = fullfile(parent, [name '.mat']);
elseif ~strcmpi(ext, '.mat')
    error('MCMC_meanage_continuous:InvalidPathway', ...
        'pathway must name a MAT file.');
end
if isempty(parent)
    parent = pwd;
    pathway = fullfile(parent, [name '.mat']);
end
if ~isfolder(parent)
    [ok, message] = mkdir(parent);
    if ~ok
        error('MCMC_meanage_continuous:OutputDirectory', '%s', message);
    end
end
paths.snapshot = pathway;
paths.checkpoint = [pathway '.checkpoint.mat'];
paths.chunks = [pathway '.chunks'];
paths.pause = [pathway '.pause'];
paths.lock = [pathway '.lock'];
if ~isfolder(paths.chunks)
    [ok, message] = mkdir(paths.chunks);
    if ~ok
        error('MCMC_meanage_continuous:ChunkDirectory', '%s', message);
    end
end
end

function cleanup = acquireLock(lockDirectory, forceUnlock)
if isfolder(lockDirectory)
    if ~forceUnlock
        error('MCMC_meanage_continuous:OutputLocked', ...
            ['Another writer may be using this output. If a previous MATLAB ' ...
             'process terminated, inspect the output and set mc.force_unlock=true: %s'], ...
            lockDirectory);
    end
    rmdir(lockDirectory);
end
[ok, message] = mkdir(lockDirectory);
if ~ok
    error('MCMC_meanage_continuous:LockFailed', '%s', message);
end
cleanup = onCleanup(@() releaseLock(lockDirectory));
end

function releaseLock(lockDirectory)
if isfolder(lockDirectory)
    rmdir(lockDirectory);
end
end

function identity = buildIdentity(mc, layout, lh_fun, settings)
identity = struct();
identity.sampler_id = 'meanage_continuous_fixed_rw_v1';
identity.proposal_kind = 'independent_uniform_original_scale';
identity.likelihood_name = func2str(lh_fun);

identity.crc_stage_reporting = 'asrstage_v1';
identity.parameter_count = layout.parameter_count;
identity.block_names = layout.block_names;
identity.active_mask = layout.active_mask;
identity.proposal_block_index = [repmat(find(layout.nh_mask),1,4), find(layout.fit_mask)];
identity.block_M = arrayfun(@(x)x.M, layout.blocks, 'UniformOutput', false);
identity.block_F = arrayfun(@(x)x.F, layout.blocks, 'UniformOutput', false);
identity.scenario_id = char(string(controlValue(mc, ...
    'continuous_scenario_id', 'unspecified')));
identity.data_id = char(string(controlValue(mc, ...
    'continuous_data_id', 'unspecified')));
identity.schedule = '4_NH_sweeps_then_1_FIT_sweep';
identity.step_definition = 'full_width_times_uniform_minus_half_to_half';
identity.tuning_block_mask = settings.tuning_block_mask;
identity.fixed_settings = selectedScenarioSettings(mc);
identity.chunk_size = settings.chunk_size;
identity.tail_size = settings.tail_size;
end

function selected = selectedScenarioSettings(mc)
names = {'contractor_list', 'No_CS_pathway_pro', 'fnc_phi', ...
    'SL_fix_target', 'SL_prior_mode', 'knot_incidence_rate', ...
    'knot_age_effect_disease_progression', 'copy_log_mode', ...
    'disable_progression_age_effect', 'inc_control', 'adj_c', ...
    'faster_odds', 'cdf_constraint', 'SL_constraint', 'Prevalence_constraint', ...
    'fnc_fit_threshold', 'fnc_stage_boost', 'rw_age_group_weights', 'rw_age_group_edges', ...
    'olrw_age_sex_weights', 'olrw_age_group_edges', 'olrw_target_lab', ...
    'knot_omega', 'num_age_eff_mean', 'num_mean_v', 'num_std_v', 'num_pZero'};
selected = struct();
for k = 1:numel(names)
    if isfield(mc, names{k})
        selected.(names{k}) = mc.(names{k});
    end
end
end

function runtime = initializeRuntime(initialM, initialF, mc, ...
        observationsM, observationsF, inmM, inmF, layout, settings, ...
        identity, lh_fun)
state = struct();
state.current = syncMCMCParameters([initialM, initialF], layout.parameter_count);
if any(~isfinite(state.current(layout.independent_indices))) || ...
        any(state.current(layout.independent_indices) < 0 | ...
            state.current(layout.independent_indices) > 1)
    error('MCMC_meanage_continuous:InvalidInitialState', ...
        'Every independently sampled initial parameter must lie in [0,1].');
end
state.mc = mc;
state.mc.copy_use_prev_cache = true;
if isfield(state.mc, 'copy_prev_cache')
    state.mc = rmfield(state.mc, 'copy_prev_cache');
end
[initialOutput, state.logL] = evaluateLikelihood(state.current, ...
    layout.parameter_count, state.mc, observationsM, observationsF, ...
    inmM, inmF, lh_fun);
if ~isfinite(state.logL)
    error('MCMC_meanage_continuous:InvalidInitialTarget', ...
        'The initial log-likelihood must be finite.');
end
state.mc = updateCache(state.mc, initialOutput);

runtime = struct();
runtime.version = 1;
runtime.identity = identity;
runtime.stage = 'tuning';
runtime.status = 'tuning_running';
runtime.state = state;
runtime.step_size = settings.initial_step_size;
runtime.total_production_samples = 0;
runtime.total_accept = zeros(1, numel(layout.blocks));
runtime.total_attempt = zeros(1, numel(layout.blocks));
runtime.report_accept = zeros(1, numel(layout.blocks));
runtime.report_attempt = zeros(1, numel(layout.blocks));
runtime.acceptance_history = emptyAcceptanceHistory();
runtime.diagnostic_history = emptyDiagnosticHistory();
runtime.tuning = initializeTuning(settings, layout);
runtime.pending = emptyBuffer(settings.chunk_size, layout.parameter_count, ...
    numel(layout.blocks), layout.proposals_per_macro);
runtime.tail = emptyBuffer(settings.tail_size, layout.parameter_count, ...
    numel(layout.blocks), layout.proposals_per_macro);
runtime.rng_state = rng;
if settings.skip_tuning
    runtime.tuning.completed = true;
    runtime.tuning.passed = true;
    runtime.tuning.skipped = true;
    runtime.tuning.reason = 'compatible_fixed_step_supplied';
    runtime.stage = 'production';
    runtime.status = 'production_running';
end
end

function tuning = initializeTuning(settings, layout)
maxCycles = settings.tuning_max_cycles;
tuning = struct();
tuning.completed = false;
tuning.passed = false;
tuning.skipped = false;
tuning.reason = '';
tuning.cycles = 0;
tuning.consecutive_passing_windows = 0;
tuning.anss1 = zeros(maxCycles, layout.parameter_count);
tuning.anss2 = zeros(maxCycles, layout.parameter_count);
tuning.likelihood_record = zeros(maxCycles, 1);
tuning.block_accept_record = zeros(maxCycles, numel(layout.blocks), 'uint8');
tuning.block_attempt_record = zeros(maxCycles, numel(layout.blocks), 'uint8');
tuning.proposal_accept_record = false(maxCycles, layout.proposals_per_macro);
nWindows = settings.tuning_max_cycles/settings.tuning_window;
tuning.window_end = zeros(nWindows, 1);
tuning.window_acceptance = nan(nWindows, numel(layout.blocks));
tuning.window_step_size = nan(nWindows, numel(layout.blocks));
tuning.window_count = 0;
end

function runtime = runTuning(runtime, observationsM, observationsF, ...
        inmM, inmF, layout, settings, identity, lh_fun, paths)
if runtime.tuning.completed
    runtime.stage = 'production';
    runtime.status = 'production_running';
    return
end

fprintf('Finite step-size calibration: cycles %d-%d, window %d.\n', ...
    settings.tuning_min_cycles, settings.tuning_max_cycles, ...
    settings.tuning_window);
lastCommit = tic;
for t = runtime.tuning.cycles + 1:settings.tuning_max_cycles
    [runtime.state, accepted, attempted, decisions] = runMacroCycle( ...
        runtime.state, layout, runtime.step_size, lh_fun, ...
        observationsM, observationsF, inmM, inmF);
    runtime.tuning.cycles = t;
    runtime.tuning.anss1(t,:) = runtime.state.current(1:layout.parameter_count);
    runtime.tuning.anss2(t,:) = runtime.state.current(layout.parameter_count+1:end);
    runtime.tuning.likelihood_record(t) = runtime.state.logL;
    runtime.tuning.block_accept_record(t,:) = uint8(accepted);
    runtime.tuning.block_attempt_record(t,:) = uint8(attempted);
    runtime.tuning.proposal_accept_record(t,:) = decisions;

    if isfile(paths.pause)
        runtime.status = 'tuning_paused';
        runtime.rng_state = rng;
        commitState(runtime, paths, settings);
        delete(paths.pause);
        fprintf('Paused safely after calibration cycle %d.\n', t);
        return
    end

    if toc(lastCommit) >= settings.save_interval_sec
        commitState(runtime, paths, settings);
        lastCommit = tic;
    end
    if mod(t, settings.tuning_window) ~= 0
        continue
    end
    first = t-settings.tuning_window+1;
    accepts = sum(double(runtime.tuning.block_accept_record(first:t,:)), 1);
    attempts = sum(double(runtime.tuning.block_attempt_record(first:t,:)), 1);
    rates = nan(1, numel(layout.blocks));
    rates(layout.active_mask) = accepts(layout.active_mask)./ ...
        attempts(layout.active_mask);
    w = runtime.tuning.window_count + 1;
    runtime.tuning.window_count = w;
    runtime.tuning.window_end(w) = t;
    runtime.tuning.window_acceptance(w,:) = rates;
    runtime.tuning.window_step_size(w,:) = runtime.step_size;
    tuned = settings.tuning_block_mask;
    passing = all(rates(tuned) >= settings.accept_min & ...
        rates(tuned) <= settings.accept_max);
    if t >= settings.tuning_min_cycles && passing
        runtime.tuning.consecutive_passing_windows = ...
            runtime.tuning.consecutive_passing_windows + 1;
    else
        runtime.tuning.consecutive_passing_windows = 0;
    end
    fprintf('  calibration cycle %d: tuned acceptance range %.3f-%.3f\n', ...
        t, min(rates(tuned)), max(rates(tuned)));

    if runtime.tuning.consecutive_passing_windows >= 2
        runtime.tuning.completed = true;
        runtime.tuning.passed = true;
        runtime.tuning.reason = 'two_consecutive_windows_passed';
        break
    end
    if t < settings.tuning_max_cycles
        high = tuned & rates > settings.accept_max;
        low = tuned & rates < settings.accept_min;
        runtime.step_size(high) = runtime.step_size(high)*settings.scale_up;
        runtime.step_size(low) = runtime.step_size(low)*settings.scale_down;
    end
    runtime.rng_state = rng;
    commitState(runtime, paths, settings);
end

if ~runtime.tuning.completed
    runtime.tuning.completed = true;
    runtime.tuning.passed = false;
    runtime.tuning.reason = 'maximum_cycles_reached';
    warning('MCMC_meanage_continuous:TuningNotPassed', ...
        ['Calibration reached %d cycles. The last tested step sizes will ' ...
         'be frozen and production will continue.'], settings.tuning_max_cycles);
end
runtime.tuning = trimTuning(runtime.tuning);
runtime.stage = 'production';
runtime.status = 'production_running';
runtime.rng_state = rng;
runtime.identity = identity;
commitState(runtime, paths, settings);
end

function tuning = trimTuning(tuning)
n = tuning.cycles;
w = tuning.window_count;
tuning.anss1 = tuning.anss1(1:n,:);
tuning.anss2 = tuning.anss2(1:n,:);
tuning.likelihood_record = tuning.likelihood_record(1:n,:);
tuning.block_accept_record = tuning.block_accept_record(1:n,:);
tuning.block_attempt_record = tuning.block_attempt_record(1:n,:);
tuning.proposal_accept_record = tuning.proposal_accept_record(1:n,:);
tuning.window_end = tuning.window_end(1:w,:);
tuning.window_acceptance = tuning.window_acceptance(1:w,:);
tuning.window_step_size = tuning.window_step_size(1:w,:);
end

function [state, accepted, attempted, decisions] = runMacroCycle(state, layout, ...
        stepSize, lh_fun, observationsM, observationsF, inmM, inmF)
accepted = zeros(1, numel(layout.blocks));
attempted = zeros(1, numel(layout.blocks));
decisions = false(1, layout.proposals_per_macro);
proposalIndex = 0;
for sweep = 1:4
    for b = find(layout.nh_mask)
        attempted(b) = attempted(b) + 1;
        [state, didAccept] = attemptBlock(state, layout, b, stepSize(b), ...
            lh_fun, observationsM, observationsF, inmM, inmF);
        accepted(b) = accepted(b) + didAccept;
        proposalIndex = proposalIndex + 1;
        decisions(proposalIndex) = didAccept;
    end
end
for b = find(layout.fit_mask)
    attempted(b) = attempted(b) + 1;
    [state, didAccept] = attemptBlock(state, layout, b, stepSize(b), ...
        lh_fun, observationsM, observationsF, inmM, inmF);
    accepted(b) = accepted(b) + didAccept;
    proposalIndex = proposalIndex + 1;
    decisions(proposalIndex) = didAccept;
end
end

function [state, accepted] = attemptBlock(state, layout, blockIndex, width, ...
        lh_fun, observationsM, observationsF, inmM, inmF)
nP = layout.parameter_count;
block = layout.blocks(blockIndex);
indices = [block.M, nP + block.F];
proposal = state.current;
jump = width*(rand(1, numel(indices))-0.5);
proposal(indices) = reflectUnitInterval(state.current(indices)+jump);
proposal = syncMCMCParameters(proposal, nP);
[candidateOutput, candidateLogL] = evaluateLikelihood(proposal, nP, ...
    state.mc, observationsM, observationsF, inmM, inmF, lh_fun);
accepted = false;
if isfinite(candidateLogL) && ...
        log(rand) < min(0, candidateLogL-state.logL)
    state.current = proposal;
    state.logL = candidateLogL;
    state.mc = updateCache(state.mc, candidateOutput);
    accepted = true;
end
end

function values = reflectUnitInterval(values)
values = mod(values, 2);
negative = values < 0;
values(negative) = values(negative) + 2;
values = 1-abs(1-values);
end

function [output, logL] = evaluateLikelihood(parameters, nP, mc, ...
        observationsM, observationsF, inmM, inmF, lh_fun)
output = lh_fun(parameters(1:nP), parameters(nP+1:end), ...
    observationsM, observationsF, mc, inmM, inmF);
if ~isstruct(output) || ~isfield(output, 'LogL') || ...
        ~isscalar(output.LogL) || ~isreal(output.LogL)
    error('MCMC_meanage_continuous:InvalidLikelihoodOutput', ...
        'The likelihood must return a struct containing a real scalar LogL.');
end
logL = double(output.LogL);
end

function mc = updateCache(mc, output)
if isfield(output, 'copy_prev_cache')
    mc.copy_prev_cache = output.copy_prev_cache;
end
end

function runtime = appendProduction(runtime, accepted, attempted, decisions, elapsed, ...
        settings, nP)
p = runtime.pending.count + 1;
runtime.pending.count = p;
runtime.pending.anss1(p,:) = runtime.state.current(1:nP);
runtime.pending.anss2(p,:) = runtime.state.current(nP+1:end);
runtime.pending.likelihood_record(p) = runtime.state.logL;
runtime.pending.time_per_sample(p) = elapsed;
runtime.pending.block_accept_record(p,:) = uint8(accepted);
runtime.pending.block_attempt_record(p,:) = uint8(attempted);
runtime.pending.proposal_accept_record(p,:) = decisions;
runtime.total_production_samples = runtime.total_production_samples + 1;
runtime.pending.global_sample_index(p) = runtime.total_production_samples;
runtime.total_accept = runtime.total_accept + accepted;
runtime.total_attempt = runtime.total_attempt + attempted;
runtime.report_accept = runtime.report_accept + accepted;
runtime.report_attempt = runtime.report_attempt + attempted;
runtime.rng_state = rng;
if p > settings.chunk_size
    error('MCMC_meanage_continuous:PendingOverflow', ...
        'Internal pending buffer exceeded chunk_size.');
end
end

function runtime = reportAcceptance(runtime, layout)
rates = nan(1, numel(layout.blocks));
rates(layout.active_mask) = runtime.report_accept(layout.active_mask)./ ...
    runtime.report_attempt(layout.active_mask);
record = struct('sample', runtime.total_production_samples, 'rates', rates);
runtime.acceptance_history(end+1) = record;
fprintf('  production sample %d, LogL %.4f, acceptance %.3f-%.3f\n', ...
    runtime.total_production_samples, runtime.state.logL, ...
    min(rates(layout.active_mask)), max(rates(layout.active_mask)));
runtime.report_accept(:) = 0;
runtime.report_attempt(:) = 0;
end

function runtime = updateDiagnostics(runtime, layout, settings)
[chainM, chainF, ~, ~, ~, ~, indices] = currentTail(runtime, settings);
fullChain = [chainM, chainF];
freeChain = fullChain(:, layout.independent_indices);
ess = computeESS(freeChain);

fraction = mean(isfinite(ess) & ess >= settings.ess_target);
minEss = min(ess);
quantiles = prctile(freeChain, [2.5, 50, 97.5], 1);
record = struct();
record.sample = runtime.total_production_samples;
record.window_start = indices(1);
record.window_end = indices(end);
record.ess = ess;
record.min_ess = minEss;
record.fraction_at_target = fraction;
record.ess_milestone = fraction >= settings.ess_stop_pct;
record.multichain_passed = false;
record.usable_for_inference = false;
if isfield(runtime.state.mc, 'convergence_paths') && ...
        numel(runtime.state.mc.convergence_paths) >= 2
    chainReport = checkMCMCConvergence(runtime.state.mc.convergence_paths);
    record.multichain_passed = chainReport.passed;
    record.usable_for_inference = chainReport.passed;
end
record.mean = mean(freeChain, 1);
record.median = quantiles(2,:);
record.ci95 = quantiles([1, 3],:);
runtime.diagnostic_history(end+1) = record;
if numel(runtime.diagnostic_history) > settings.max_diagnostic_records
    runtime.diagnostic_history = runtime.diagnostic_history( ...
        end-settings.max_diagnostic_records+1:end);
end
fprintf('    diagnostics: min ESS %.1f; %.1f%% at ESS >= %.0f%s\n', ...
    minEss, 100*fraction, settings.ess_target, ...
    milestoneText(record.ess_milestone));
end

function text = milestoneText(reached)
if reached
    text = ' [single-chain ESS milestone only; confirm independent chains]';
else
    text = '';
end
end

function ess = computeESS(chain)
[n, p] = size(chain);
ess = zeros(1, p);
if n < 3
    return
end
nFft = 2^nextpow2(2*n-1);
for k = 1:p
    values = chain(:,k);
    if all(values==values(1)), continue; end
    values = values-mean(values);
    variance = sum(values.^2)/n;
    if ~isfinite(variance) || variance <= 0
        continue
    end
    spectrum = fft(values, nFft);
    autocovariance = ifft(abs(spectrum).^2, 'symmetric');
    autocovariance = autocovariance(1:n)./(n-(0:n-1)');
    rho = autocovariance/variance;
    pairSum = rho(2:2:end-1)+rho(3:2:end);
    stop = find(pairSum <= 0 | ~isfinite(pairSum), 1, 'first');
    if isempty(stop)
        kept = pairSum;
    else
        kept = pairSum(1:stop-1);
    end
    tau = max(1, 1+2*sum(kept));
    ess(k) = min(n, n/tau);
end
end

function runtime = flushPending(runtime, paths, settings)
n = runtime.pending.count;
if n == 0
    return
end
chunk = trimBuffer(runtime.pending, n);
chunk.global_start = chunk.global_sample_index(1);
chunk.global_end = chunk.global_sample_index(end);
chunkFile = fullfile(paths.chunks, sprintf('chunk_%012d_%012d.mat', ...
    chunk.global_start, chunk.global_end));
if isfile(chunkFile)
    old = load(chunkFile, 'chunk');
    if ~isfield(old, 'chunk') || old.chunk.global_start ~= chunk.global_start || ...
            old.chunk.global_end ~= chunk.global_end || ...
            ~isequaln(old.chunk.anss1, chunk.anss1) || ...
            ~isequaln(old.chunk.anss2, chunk.anss2)
        error('MCMC_meanage_continuous:ChunkConflict', ...
            'Existing chunk conflicts with the pending samples: %s', chunkFile);
    end
else
    atomicSave(chunkFile, 'chunk', chunk);
end
runtime.tail = appendTail(runtime.tail, chunk, settings.tail_size);
runtime.pending = emptyBuffer(settings.chunk_size, ...
    runtime.identity.parameter_count, numel(runtime.step_size), ...
    numel(runtime.identity.proposal_block_index));
end

function tail = appendTail(tail, chunk, tailSize)
old = trimBuffer(tail, tail.count);
combined = struct();
fields = {'anss1', 'anss2', 'likelihood_record', 'time_per_sample', ...
    'block_accept_record', 'block_attempt_record', 'global_sample_index', ...
    'proposal_accept_record'};
for k = 1:numel(fields)
    combined.(fields{k}) = [old.(fields{k}); chunk.(fields{k})];
end
keepStart = max(1, size(combined.anss1,1)-tailSize+1);
tail = emptyBuffer(tailSize, size(combined.anss1,2), ...
    size(combined.block_accept_record,2), size(combined.proposal_accept_record,2));
tail.count = size(combined.anss1,1)-keepStart+1;
for k = 1:numel(fields)
    values = combined.(fields{k})(keepStart:end,:);
    tail.(fields{k})(1:tail.count,:) = values;
end
end

function buffer = emptyBuffer(capacity, nP, numBlocks, numProposals)
buffer = struct();
buffer.count = 0;
buffer.anss1 = zeros(capacity, nP);
buffer.anss2 = zeros(capacity, nP);
buffer.likelihood_record = zeros(capacity, 1);
buffer.time_per_sample = zeros(capacity, 1);
buffer.block_accept_record = zeros(capacity, numBlocks, 'uint8');
buffer.block_attempt_record = zeros(capacity, numBlocks, 'uint8');
buffer.global_sample_index = zeros(capacity, 1);
buffer.proposal_accept_record = false(capacity, numProposals);
end

function trimmed = trimBuffer(buffer, n)
fields = {'anss1', 'anss2', 'likelihood_record', 'time_per_sample', ...
    'block_accept_record', 'block_attempt_record', 'global_sample_index', ...
    'proposal_accept_record'};
trimmed = struct();
trimmed.count = n;
for k = 1:numel(fields)
    trimmed.(fields{k}) = buffer.(fields{k})(1:n,:);
end
end

function [M, F, logL, elapsed, accepted, attempted, indices, decisions] = ...
        currentTail(runtime, settings)
tail = trimBuffer(runtime.tail, runtime.tail.count);
pending = trimBuffer(runtime.pending, runtime.pending.count);
fields = {'anss1', 'anss2', 'likelihood_record', 'time_per_sample', ...
    'block_accept_record', 'block_attempt_record', 'global_sample_index', ...
    'proposal_accept_record'};
values = cell(1, numel(fields));
for k = 1:numel(fields)
    values{k} = [tail.(fields{k}); pending.(fields{k})];
end
start = max(1, size(values{1},1)-settings.tail_size+1);
for k = 1:numel(values)
    values{k} = values{k}(start:end,:);
end
[M, F, logL, elapsed, accepted, attempted, indices, decisions] = values{:};
end

function commitState(runtime, paths, settings)
runtime.rng_state = rng;
atomicSave(paths.checkpoint, 'runtime', runtime);
answer = buildAnswer(runtime, paths, settings);
atomicSave(paths.snapshot, 'answer', answer);
end

function answer = buildAnswer(runtime, paths, settings)
[M, F, logL, elapsed, accepted, attempted, indices, decisions] = ...
    currentTail(runtime, settings);
answer = struct();
answer.sampler_id = runtime.identity.sampler_id;
answer.proposal_kind = runtime.identity.proposal_kind;
answer.parameter_count = runtime.identity.parameter_count;
answer.likelihood_name = runtime.identity.likelihood_name;
answer.scenario_id = runtime.identity.scenario_id;
answer.data_id = runtime.identity.data_id;
answer.block_names = runtime.identity.block_names;
answer.block_active = runtime.identity.active_mask;
answer.block_M = runtime.identity.block_M;
answer.block_F = runtime.identity.block_F;
answer.anss1 = M;
answer.anss2 = F;
answer.likelihood_record = logL;
answer.time_per_sample = elapsed;
answer.block_accept_record = accepted;
answer.block_attempt_record = attempted;
answer.proposal_accept_record = decisions;
answer.proposal_block_index = runtime.identity.proposal_block_index;
answer.usable_for_inference = ~isempty(runtime.diagnostic_history) && ...
    runtime.diagnostic_history(end).usable_for_inference;
answer.global_sample_index = indices;
answer.total_production_samples = runtime.total_production_samples;
answer.tail_start_sample = firstOrZero(indices);
answer.parmStepSize_record = runtime.step_size;
answer.production_acceptance = safeRates( ...
    runtime.total_accept, runtime.total_attempt);
answer.paccept_temp = acceptanceHistoryMatrix(runtime.acceptance_history, numel(runtime.step_size));
answer.tuning = runtime.tuning;
answer.diagnostics = runtime.diagnostic_history;
answer.status = runtime.status;
answer.stage = runtime.stage;
answer.production_complete = strcmp(runtime.status, 'production_complete');
answer.final_LogL = runtime.state.logL;
answer.checkpoint_file = paths.checkpoint;
answer.chunk_directory = paths.chunks;
answer.settings = settings;
answer.identity = runtime.identity;
end

function value = firstOrZero(values)
if isempty(values)
    value = 0;
else
    value = values(1);
end
end

function rates = safeRates(accepts, attempts)
rates = nan(size(accepts));
used = attempts > 0;
rates(used) = accepts(used)./attempts(used);
end

function matrix = acceptanceHistoryMatrix(history, numBlocks)
if isempty(history)
    matrix = zeros(0, numBlocks+1);
    return
end
matrix = zeros(numel(history), numBlocks+1);
for k = 1:numel(history)
    matrix(k,:) = [history(k).sample, history(k).rates];
end
end

function history = emptyAcceptanceHistory()
history = struct('sample', {}, 'rates', {});
end

function history = emptyDiagnosticHistory()
history = struct('sample', {}, 'window_start', {}, 'window_end', {}, ...
    'ess', {}, 'min_ess', {}, 'fraction_at_target', {}, ...
    'ess_milestone', {}, 'multichain_passed', {}, 'usable_for_inference', {}, 'mean', {}, 'median', {}, 'ci95', {});
end

function runtime = loadCheckpoint(checkpointFile, identity, settings)
if ~isfile(checkpointFile)
    error('MCMC_meanage_continuous:MissingCheckpoint', ...
        'Resume requested but checkpoint does not exist: %s', checkpointFile);
end
saved = load(checkpointFile, 'runtime');
if ~isfield(saved, 'runtime') || ~isstruct(saved.runtime)
    error('MCMC_meanage_continuous:InvalidCheckpoint', ...
        'Checkpoint does not contain a runtime struct.');
end
runtime = saved.runtime;
runtime.identity = normalizeStageNames(runtime.identity);
if isfield(runtime.state.mc, 'copy_prev_cache') && ...
        ~isfield(runtime.state.mc.copy_prev_cache, 'state_names')
    runtime.state.mc = rmfield(runtime.state.mc, 'copy_prev_cache');
end
runtime.state.mc = normalizeStageNames(runtime.state.mc);
if isfield(runtime, 'identity') && isfield(runtime.identity, 'likelihood_name') && ...
        strcmp(runtime.identity.likelihood_name, 'loglikelihood_twosample_nh_meanage_base') && ...
        strcmp(identity.likelihood_name, 'logLikelihood_CRC')
    runtime.identity.likelihood_name = identity.likelihood_name;
end
if ~isfield(runtime, 'identity') || ~isequaln(runtime.identity, identity)
    error('MCMC_meanage_continuous:IdentityMismatch', ...
        ['Checkpoint identity does not match the requested likelihood, ' ...
         'scenario, data, parameter layout, or fixed settings.']);
end
if size(runtime.pending.anss1, 1) ~= settings.chunk_size || ...
        size(runtime.tail.anss1, 1) ~= settings.tail_size
    error('MCMC_meanage_continuous:StorageSettingMismatch', ...
        'chunk_size and tail_size must match the checkpoint.');
end
end

function atomicSave(path, variableName, value)
temporary = [path '.tmp'];
payload = struct();
payload.(variableName) = value;
save(temporary, '-struct', 'payload', '-v7.3');
[ok, message] = movefile(temporary, path, 'f');
if ~ok
    error('MCMC_meanage_continuous:AtomicSaveFailed', ...
        'Could not replace %s: %s', path, message);
end
end

function value = controlValue(mc, fieldName, defaultValue)
if isfield(mc, fieldName) && ~isempty(mc.(fieldName))
    value = mc.(fieldName);
else
    value = defaultValue;
end
end

function value = logicalControl(mc, fieldName, defaultValue)
value = controlValue(mc, fieldName, defaultValue);
if ~isscalar(value) || ~(islogical(value) || isnumeric(value))
    error('MCMC_meanage_continuous:InvalidLogicalControl', ...
        'mc.%s must be a logical scalar.', fieldName);
end
value = logical(value);
end
