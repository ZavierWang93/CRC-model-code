function value = normalizeStageNames(value)

legacy = {'HP','SL','inc_SL','num_inc_SL','inc_SL_20_45','inc_SL_70', ...
    'CRC_HSL','incidence_CRC_HSL_1y','age_specific_prop_CRC_HSL','overall_prop_CRC_HSL', ...
    'HSL_constraint','HSL_fix_target','HSL_prior_mode'};
current = {'LSL','HSL','inc_LSL','num_inc_LSL','inc_LSL_20_45','inc_LSL_70', ...
    'CRC_SL','incidence_CRC_SL_1y','age_specific_prop_CRC_SL','overall_prop_CRC_SL', ...
    'SL_constraint','SL_fix_target','SL_prior_mode'};

if isstruct(value)
    names = fieldnames(value);
    for k = 1:numel(names)
        oldName = names{k};
        match = find(strcmp(oldName, legacy), 1);
        newName = oldName;
        if ~isempty(match), newName = current{match}; end
        for row = 1:numel(value)
            value(row).(newName) = normalizeStageNames(value(row).(oldName));
        end
        if ~strcmp(oldName, newName), value = rmfield(value, oldName); end
    end
elseif iscell(value)
    for k = 1:numel(value), value{k} = normalizeStageNames(value{k}); end
elseif ischar(value)
    match = find(strcmp(value, legacy), 1);
    if ~isempty(match), value = current{match}; end
elseif isstring(value)
    for k = 1:numel(legacy), value(value == legacy{k}) = string(current{k}); end
end
end
