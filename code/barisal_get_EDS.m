function [EDS, eds_i] = barisal_get_EDS(EDS,entity_file,hazard_file)

for eds_i =1:length(EDS)
    if strcmp(EDS(eds_i).assets.filename,entity_file) ...
            && strcmp(EDS(eds_i).hazard.filename,hazard_file)
        EDS = EDS(eds_i);
        return
    end
end
fprintf('ERROR: EDS not found\n')
EDS = []; eds_i = [];
return