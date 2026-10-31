function layout = mcmcBlockLayout(nP, mc)

names = {'prevalence','asr','inc_NA','inc_LSL_20_45','inc_LSL_70', ...
    'exit_rate','age_effect','prob_diagnosed','mean_v_lab1', ...
    'std_v_lab1','pZero_lab1','q0_lab1','qNA_lab1','omega', ...
    'age_eff_mean','mean_v_lab2','std_v_lab2','pZero_lab2','q0_lab2','qNA_lab2'};
if ~ismember(nP,[64 99 105 106 113])
    error('MeanageLayout:ParameterCount','Unsupported meanage layout.');
end
nFit = 6;
if nP == 113, nFit = 7; end
if nP == 99
    nh = {1:5,6:14,15,16,[],17:24,[],25:27}; nhEnd = 27;
else
    nh = {1:5,6:14,15:17,18:19,20,21:28,29:30,31:33}; nhEnd = 33;
end
lab1 = nhEnd + 1;
lab2 = lab1 + 4*nFit + 7;
indices = [nh, {lab1+(0:nFit-1), lab1+nFit+(0:nFit-1), ...
    lab1+2*nFit+(0:nFit-1),lab1+3*nFit,lab1+3*nFit+1, ...
    lab1+3*nFit+(2:6),lab1+3*nFit+7, ...
    lab2+(0:nFit-1),lab2+nFit+(0:nFit-1), ...
    lab2+2*nFit+(0:nFit-1),lab2+3*nFit,lab2+3*nFit+1}];
active = ~cellfun(@isempty, indices);
if nP == 64
    if numel(mc.contractor_list) ~= 1 || ~ismember(mc.contractor_list,[1 2])
        error('MeanageLayout:Lab','The 64-value model needs one selected lab.');
    end
    active(16:20) = false;
    indices(16:20) = {[],[],[],[],[]};
    if mc.contractor_list == 2
        names(9:13) = {'mean_v_lab2','std_v_lab2','pZero_lab2','q0_lab2','qNA_lab2'};
    end
end
if nP == 106
    names{21} = 'fraction_B';
    indices{21} = 106;
    active(21) = true;
end

if isfield(mc,'num_mean_v')
    actual = [mc.num_prevalence mc.num_ascertainment_ratio mc.num_inc_NA ...
        mc.num_inc_LSL mc.num_exit_rate mc.num_age_effect_on_disease_progression ...
        mc.num_prob_diagnosed mc.num_mean_v mc.num_std_v mc.num_pZero ...
        mc.num_age_eff_mean];
    expected = [6 9 3 3 8 2 3 nFit nFit nFit nFit];
    if nP == 99, expected(3:6) = [1 1 8 0]; end
    if ~isequal(actual,expected)
        error('MeanageLayout:Counts','Model counts do not match the selected layout.');
    end
end
blocks = repmat(struct('name','','M',[],'F',[],'active',false,'group',''),1,numel(names));
nhMask = false(size(active)); nhMask(1:8) = active(1:8);
if nP == 106, nhMask(21) = true; end
allIndices = [];
for k = 1:numel(names)
    blocks(k).name = names{k};
    blocks(k).M = indices{k};
    blocks(k).F = indices{k};
    if k == 14 || k == 21, blocks(k).F = []; end
    blocks(k).active = active(k);
    blocks(k).group = 'FIT';
    if nhMask(k), blocks(k).group = 'NH'; end
    if active(k)
        allIndices = [allIndices blocks(k).M nP+blocks(k).F];
    end
end
assert(numel(allIndices)==numel(unique(allIndices)) && all(allIndices<=2*nP), ...
    'Independent parameter indices must be unique and in range.');
layout = struct('parameter_count',nP,'blocks',blocks,'block_names',{names}, ...
    'active_mask',active,'nh_mask',nhMask,'fit_mask',active & ~nhMask, ...
    'independent_indices',sort(allIndices),'independent_parameter_count',numel(allIndices), ...
    'proposals_per_macro',4*nnz(nhMask)+nnz(active & ~nhMask));
end
