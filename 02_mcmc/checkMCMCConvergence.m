function report = checkMCMCConvergence(pathways)

pathways = cellstr(string(pathways));
if numel(pathways)<2 || numel(unique(pathways))~=numel(pathways)
    error('MeanageConvergence:Paths','Supply >=2 different independent-chain paths.');
end
chains = cell(size(pathways)); identities = cell(size(pathways));
sampleIndices = cell(size(pathways));
for c=1:numel(pathways)
    if ~isfile(pathways{c})
        report=struct('passed',false,'reason','Waiting for all chain snapshots.');
        fprintf('  [Convergence] %s\n',report.reason); return;
    end
    saved=load(pathways{c},'answer');
    if ~isfield(saved,'answer') || ~isfield(saved.answer,'identity')
        error('MeanageConvergence:Sampler','Use snapshots from the new meanage sampler.');
    end
    a=saved.answer; identities{c}=normalizeStageNames(a.identity);
    sampleIndices{c}=a.global_sample_index;
    indices=[];
    for b=find(a.block_active)
        indices=[indices a.block_M{b} a.parameter_count+a.block_F{b}];
    end
    values=[a.anss1 a.anss2];
    chains{c}=values(:,sort(indices));
    if c>1 && isequaln(chains{c},chains{1})
        error('MeanageConvergence:Duplicate','Copied or identical chains are not independent runs.');
    end

    if c>1 && ~isequaln(identities{c},identities{1})
        error('MeanageConvergence:Identity','These chains do not have the same model/data/settings.');
    end
end
n=min(cellfun(@(x)size(x,1),chains));
if n<1000
    report=struct('passed',false,'reason','Need >=1000 saved production draws per chain.');
    fprintf('  [Convergence] %s\n',report.reason); return;
end
draws=zeros(n,numel(chains),size(chains{1},2));
for c=1:numel(chains), draws(:,c,:)=reshape(chains{c}(end-n+1:end,:),n,1,[]); end
report=chainDiagnostics(draws);
report.paths=pathways;
report.window_start=zeros(size(pathways));
report.window_end=zeros(size(pathways));
for c=1:numel(pathways)
    indices=sampleIndices{c};
    report.window_start(c)=indices(end-n+1);
    report.window_end(c)=indices(end);
end
if report.passed
    report.reason='Independent-chain diagnostic thresholds passed; inspect traces/initialization before reporting.';
else
    report.reason='Not yet usable: R-hat or bulk/tail ESS thresholds failed.';
end
fprintf('  [Convergence] %d chains: max rank R-hat %.4f, min bulk ESS %.1f, min tail ESS %.1f.\n', ...
    report.num_chains,report.max_rhat,report.min_bulk_ess,report.min_tail_ess);
fprintf('  [Convergence] %s Sampling continues.\n',report.reason);
end
