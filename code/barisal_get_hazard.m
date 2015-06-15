function [hazard, h_i] = barisal_get_hazard(reference_year,cc,peril_ID,hazard_files)
for h_i = 1:length(hazard_files)
    if isempty(strfind(hazard_files{h_i},num2str(reference_year)))
        continue
    end
    if ~isempty(cc) && isempty(strfind(hazard_files{h_i},cc))
        continue
    end
    load(hazard_files{h_i})
    if strcmp(hazard.peril_ID,peril_ID)
        return;
    end
end
fprintf('ERROR: specified hazard not found\n')
hazard = []; h_i = [];
end