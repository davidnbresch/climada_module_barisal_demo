function [hazard, entity, label] = barisal_hazard_entity_load(hazard_name, cc_scenario, timehorizon)
% load hazard and entity for a given hazard, timehorizon and cc_scenario
% MODULE: 
%   barisal_demo
% NAME:
%   barisal_hazard_entity_load
% PURPOSE:
%   load hazard and entity for a given hazard, timehorizon and cc_scenario
% CALLING SEQUENCE:
%   [hazard, entity, label] = barisal_hazard_entity_load(hazard_name, cc_scenario, timehorizon)
% EXAMPLE:
%   [hazard, entity, label] = barisal_hazard_entity_load('flood_depth', 'moderate', 2030)
% INPUTS:
%   hazard_name: can either be 'flood_depth', 'flood_duration', or 'cyclone_wind'
%   cc_scenario: can either be 'no change', 'moderate', or 'extreme'
%   timehorizon: can either be 2014, 2030, or 2050
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   load hazard and entity for a given hazard, timehorizon and cc_scenario,
%   label: a structure with .hazad_name, .timehorizon, .cc_scenario and
%   .titlestr to characterize the scenario
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150429, init
%-


global climada_global
if ~climada_init_vars, return; end

% poor man's version to check arguments
if ~exist('hazard_name'     ,'var'), hazard_name = []; end
if ~exist('cc_scenario'     ,'var'), cc_scenario = []; end
if ~exist('timehorizon'     ,'var'), timehorizon = []; end

% init
hazard = [];
entity = [];

label.hazard_name  = hazard_name;
label.timehorizon = timehorizon;
label.cc_scenario  = cc_scenario;
label.titlestr     = '';


fprintf('\n-----------\nSELECTED SCENARIO: %s, %d, %s\n-----------\n', hazard_name, timehorizon, cc_scenario)

%% hazard and entity
switch hazard_name
    case 'flood_depth_monsoon'
        % hazard = climada_asci2hazard(asci_file);
        hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_depth_monsoon'];
        entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 060515_flood_depth.mat'];
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 040515_flood_depth.mat'];
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_depth.mat'];
        
    case 'flood_depth_cyclone'
        % hazard = climada_asci2hazard(asci_file);
        hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_depth_cyclone'];
        entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 060515_flood_depth.mat'];
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 040515_flood_depth.mat'];
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_depth.mat'];

    case 'flood_duration_monsoon' % hazard duration (too be prepared!)
        hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_duration_monsoon'];
        entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 060515_flood_duration.mat']
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 040515_flood_duration.mat'];
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_duration.mat'];
            
    case 'flood_duration_cyclone' % hazard duration (too be prepared!)
        hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_duration_cyclone'];
        entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 060515_flood_duration.mat'];
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 040515_flood_duration.mat'];
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_duration.mat'];

    case 'cyclone_wind'
        hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC'];
        entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Cyclones 060515_cyclone_windspeed.mat'];
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Cyclones 040515_cyclone_windspeed.mat'];
        %entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.mat'];
        %cyclone_wind   = 1; flood_depth    = 0; flood_duration = 0; 
        % wind centroids
        %centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
        %load(centroids_file)   
end

%% add extensions depending on timehorizon and cc_scenario
if timehorizon == 2014 | strcmp(cc_scenario,'no change')
    hazard_set_file = sprintf('%s_%d.mat',hazard_set_file,2014);
else
    hazard_set_file = sprintf('%s_cc_%d_%s.mat',hazard_set_file,timehorizon,cc_scenario);
end

%% check if hazard file exists
if ~exist(hazard_set_file,'file')
    fprintf('Hazard does not exist. Please check for %s, %d, %s.\n', hazard_name, timehorizon, cc_scenario)
    return
else
    load(hazard_set_file)
    indx = strfind(hazard_set_file,filesep);
    fprintf('\t - Loaded hazard: %s\n', hazard_set_file(indx(end)+1:end))
    if ~isfield(hazard,'filename')
        hazard.filename = '';
        save(hazard_set_file,'hazard')
    end
end


%% check if entity file exists   
if ~exist(entity_filename,'file')
    fprintf('Entity does not exist. Please check for %s.\n', hazard_name)
    return
else
    load(entity_filename)
    indx = strfind(entity_filename,filesep);
end
  

%% assign entity values depending on timehorizon
switch timehorizon
    case 2014
        % risk today
        %entity.assets.Value = entity_ori.assets.Value;
        entity.assets.Value = entity.assets.Value;
        label.titlestr = 'Risk today';
    case 2030
        % risk 2030
        entity.assets.Value = entity.assets.Value_2030;
        label.titlestr = 'Socio-economic 2030 (scenario 1)';
    case 2050
        % risk 2050
        entity.assets.Value = entity.assets.Value_2050;
        label.titlestr = 'Socio-economic 2050 (scenario 1)';
end
fprintf('\t - Loaded entity: %s for %d\n', entity_filename(indx(end)+1:end), timehorizon)
    
    
  
    
    
