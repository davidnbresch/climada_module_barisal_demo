function [entity, e_i] = barisal_get_entity(reference_year,peril,entity_files,eco_scen)
% barisal_get_entity
% MODULE:
%   barisal_demo
% NAME:
%   barisal_get_entity
% PURPOSE:
%   Get/load entity for barisal, based on exisiting mat-files that are
%   saved in barisal_demo/data/entities/....mat as specified in entiy_files
% CALLING SEQUENCE:
%   [entity, e_i] = barisal_get_entity(reference_year,peril,entity_files,eco_scen)
% EXAMPLE:
%   entity = barisal_get_entity(2014,'cyclones',entity_files,'scenario 1')
% INPUTS:
%   reference_year: 2014, 2030 or 2050
%   peril: 'cyclones' or 'floods'
%   entity_files: a cell with all the existing entity-mat files e.g. 
%        entity_files{1} = '...barisal_demo\data\entities\Assets_at_risk_100x100_09092015 Baseline scenario 1_cyclones_2014.mat';
%        entity_files{2} = '...barisal_demo\data\entities\Assets_at_risk_100x100_09092015 Baseline scenario 1_cyclones_2030.mat';
%   eco_scen: 'scenario 1', or 'scenario 2', or 'scenario 3'
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   entity: a climada entity structure, loaded
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150909, add documentation
% Lea Mueller, muellele@gmail.com, 20150910, load entity 2014 scenario 1 to copy damagefunctions and discount
%-

if ~exist('eco_scen','var'),eco_scen='';end
if isempty(eco_scen), eco_scen = 1; end

for e_i = 1:length(entity_files)    
    if isempty(strfind(entity_files{e_i},num2str(reference_year)))
        continue
    end
    if isempty(strfind(lower(entity_files{e_i}),lower(peril)))
        continue
    end
    if ~isempty(strfind(entity_files{e_i},sprintf('scenario %d',eco_scen)))
        %[pathstr, name, ext] = fileparts(entity_files{e_i});
        %fprintf('%s\n',name)
        load(entity_files{e_i})
        
        if ~isfield(entity,'damagefunctions')
            fprintf('Add damagefunctions and discount from scenario 1\n')
            entity_copy = entity;
            entity = barisal_get_entity(2014,peril,entity_files,1);
            entity_copy.damagefunctions = entity.damagefunctions;
            entity_copy.discount = entity.discount;
            clear entity
            entity = entity_copy;
            save(entity_files{e_i},'entity')
        end
        return
    end
    
end
fprintf('ERROR: specified entity not found\n')
entity = []; e_i = [];
end


