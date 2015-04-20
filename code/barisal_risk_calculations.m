% barisal_risk_calculations

% - tc wind hazard: use barisal_tc_hazard_prob.m to create wind hazard
% - flood hazard: read asci-file from Ruud (Witteveen+Bos), load flood hazard



% tc wind
hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_prob'];
load(hazard_set_file)

% wind centroids
centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
load(centroids_file)

% flood
% asci_file = ;
% hazard = climada_asci2hazard(asci_file);
hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL'];
load(hazard_set_file)

% load entity (asset portfolio)
entity_file = [climada_global.data_dir filesep 'entities' filesep 'Barisal_BCC_1km_100.mat'];
if exist(entity_file,'file')
    load(entity_file)
end

% read ecorys entity
entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal.xls'];
[entity,entity_save_file] = climada_entity_read(entity_filename,hazard);
% convert local coordinates to lat lon
entity.assets.X = entity.assets.lon;
entity.assets.Y = entity.assets.lat;
[entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.X, entity.assets.Y);
entity = climada_assets_encode(entity,hazard);
save(entity_save_file, 'entity')
% plot for first visual check
climada_entity_plot(entity)


% calculate damage today
annotation_name = 'Flood today';
force_re_encode = 1;
silent_mode     = 0;
EDS = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
climada_EDS_DFC(EDS)




