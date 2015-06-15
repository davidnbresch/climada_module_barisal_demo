function [entity, e_i] = barisal_get_entity(reference_year,peril,entity_files)
for e_i = 1:length(entity_files)
    if isempty(strfind(entity_files{e_i},num2str(reference_year)))
        continue
    end
    if ~isempty(strfind(lower(entity_files{e_i}),lower(peril)))
        load(entity_files{e_i})
        return
    end
end
fprintf('ERROR: specified entity not found\n')
entity = []; e_i = [];
end