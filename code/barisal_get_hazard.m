function [hazard, h_i] = barisal_get_hazard(reference_year,cc,peril_ID,hazard_files)

if iscell(peril_ID), peril_ID = char(peril_ID); end
if iscell(cc),       cc       = char(cc);       end
for h_i = 1:length(hazard_files)
    if isempty(strfind(hazard_files{h_i},num2str(reference_year)))
        continue
    end
    if ~isempty(cc) && isempty(strfind(hazard_files{h_i},cc))
        continue
    end
    if ~isempty(strfind(hazard_files{h_i},peril_ID))
        load(hazard_files{h_i})
        return;
    end
end
fprintf('ERROR: specified hazard not found\n')
hazard = []; h_i = [];
end