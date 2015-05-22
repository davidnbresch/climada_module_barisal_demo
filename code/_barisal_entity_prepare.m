%% prepare barisal entity


%% --------------------
%  FLOOD
%----------------------

flood = 1;

%% hazard
if flood == 1
    % hazard flood
    % asci_file = ;
    % hazard = climada_asci2hazard(asci_file);
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL'];
    load(hazard_set_file)
    
else    
    % hazard tc wind
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_prob'];
    load(hazard_set_file)

    % wind centroids
    centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
    load(centroids_file)
end


%% read ecorys entity Flooding
fprintf('Read entity flood\n')
entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal.xls'];
entity_filename_mat = strrep(entity_filename,'.xls', '.mat');
if exist(entity_filename_mat,'file')
    delete(entity_filename_mat)
end
[entity,entity_save_file] = climada_entity_read(entity_filename,hazard);
% convert local coordinates to lat lon
fprintf('\t- Convert to lat lon and reencode\n')
entity.assets.X = entity.assets.lon;
entity.assets.Y = entity.assets.lat;
[entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.X, entity.assets.Y);
entity = climada_assets_encode(entity,hazard);
% save(entity_save_file, 'entity')
% % plot for first visual check
% figure
% climada_entity_plot(entity,8)


%% organize damage functions (flood depth, flood duration, cyclone wind speed)

% find damagefunctions for flood duration (max intensity value is 14 days)
fprintf('\t- Organize damage functions (flood depth, flood duration, cyclone wind speed)\n')
entity_ori         = entity;
max_value_duration =  14.0;
max_value_depth    =   4.0;
max_value_windspeed= 400.0;

DamFunID       = unique(entity.damagefunctions.DamageFunID);
indx_depth     = zeros(size(DamFunID));
indx_duration  = zeros(size(DamFunID));
indx_windspeed = zeros(size(DamFunID));
for d_i = 1:length(DamFunID)
    indx = entity.damagefunctions.DamageFunID == DamFunID(d_i);
    max_value = max(entity.damagefunctions.Intensity(indx));
    
    if max_value<=max_value_depth+10^-2 & max_value>=max_value_depth-10^-2
        indx_depth(d_i) = 1;
    elseif max_value<=max_value_duration+10^-2 & max_value>=max_value_duration-10^-2
        indx_duration(d_i) = 1;
    elseif max_value<=max_value_windspeed+10^-2 & max_value>=max_value_windspeed-10^-2
        indx_windspeed(d_i) = 1;
    end
end
indx_depth     = logical(indx_depth);
indx_duration  = logical(indx_duration);
indx_windspeed = logical(indx_windspeed);

% flood depth
% find all assets that do not correspond to the specific index (damage function unit)
fprintf('\t- DamageFunctions for flood depth: %d\n', numel(DamFunID(indx_depth)))
fprintf('%d, ', DamFunID(indx_depth))
fprintf('\n')
entity           = entity_ori;
non_valid_DamFun = DamFunID(~indx_depth);
non_valid_indx   = ismember(entity.assets.DamageFunID, non_valid_DamFun);
entity.assets.Value(non_valid_indx)      = 0;
entity.assets.Value_2030(non_valid_indx) = 0;
entity.assets.Value_2050(non_valid_indx) = 0;
entity.assets.comment = 'Flood depth (m)';
entity_filename       = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_depth.mat'];
save(entity_filename, 'entity')
fprintf('\t- Save entity flood depth as \n\t%s\n\n', entity_filename)

% flood duration
% find all assets that do not correspond to the specific index (damage function unit)
fprintf('\t- DamageFunctions for flood duration: %d\n', numel(DamFunID(indx_duration)))
fprintf('%d, ', DamFunID(indx_duration))
fprintf('\n')
entity           = entity_ori;
non_valid_DamFun = DamFunID(~indx_duration);
non_valid_indx   = ismember(entity.assets.DamageFunID, non_valid_DamFun);
entity.assets.Value(non_valid_indx)      = 0;
entity.assets.Value_2030(non_valid_indx) = 0;
entity.assets.Value_2050(non_valid_indx) = 0;
entity.assets.comment = 'Flood duration (days)';
entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_duration.mat'];
save(entity_filename, 'entity')
fprintf('\t- Save entity flood duration as \n\t%s\n\n', entity_filename)




%% --------------------
%  CYCLONE WIND
%----------------------


%% hazard tc wind
hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_prob'];
load(hazard_set_file)

% wind centroids
centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
load(centroids_file)
  

%% read ecorys entity cyclones
fprintf('Read entity cyclone wind\n')
clear entity
entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.xls'];
entity_filename_mat = strrep(entity_filename,'.xls', '.mat');
if exist(entity_filename_mat,'file')
    delete(entity_filename_mat)
end
[entity,entity_save_file] = climada_entity_read(entity_filename,hazard);
% convert local coordinates to lat lon
fprintf('\t- Convert to lat lon and reencode\n')
entity.assets.X = entity.assets.lon;
entity.assets.Y = entity.assets.lat;
[entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.X, entity.assets.Y);
entity = climada_assets_encode(entity,hazard);
% save(entity_save_file, 'entity')
% plot for first visual check
figure
climada_entity_plot(entity,8)

% country_name = 'Barisal';
% check_printplot = 0;
% printname = '';
% keep_boundary = 0;
% figure
% climada_plot_entity_assets(entity,centroids,country_name,check_printplot,printname,keep_boundary);


%% organize damage functions (flood depth, flood duration, cyclone wind speed)
%  see above
%  find all assets that do not correspond to the specific index (damage function unit)
fprintf('\t- DamageFunctions for wind speed: %d\n', numel(DamFunID(indx_windspeed)))
fprintf('%d, ', DamFunID(indx_windspeed))
fprintf('\n')
non_valid_DamFun = DamFunID(~indx_windspeed);
non_valid_indx   = ismember(entity.assets.DamageFunID, non_valid_DamFun);
entity.assets.Value(non_valid_indx)      = 0;
entity.assets.Value_2030(non_valid_indx) = 0;
entity.assets.Value_2050(non_valid_indx) = 0;

% transform from kilometers per hour (kph) to m/s
entity.damagefunctions.Intensity_ori = entity.damagefunctions.Intensity;
entity.damagefunctions.Intensity     = entity.damagefunctions.Intensity_ori/3.6;
entity.damagefunctions.comment       = 'Transformed intensity to m/s from km/h';

entity.assets.comment = 'Cyclone wind speed (m/s)';
entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.mat'];
save(entity_filename, 'entity')
fprintf('\t- Save entity cyclone wind as \n\t%s\n\n', entity_filename)



%% next time only load entities
% if flood_depth ==1
%     entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_depth.mat'];
%     load(entity_filename)
%     %entity_flood = entity;
%     
% elseif flood_duration ==1
%     entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_duration.mat'];
%     load(entity_filename)
%     %entity_flood = entity;
%     
% else
%     % tc wind
%     entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.mat'];
%     load(entity_filename)
%     %entity_wind = entity;
% end



%% see damage functions for different asset categories
% 
% asset_cat = unique(entity.assets.Category(entity.assets.Value>0));
% for cat_i = 1:length(asset_cat)
%     fprintf('-----------\n-----------\nAsset category: %s \n-----------\n',asset_cat{cat_i})
%     indx = strcmp(entity.assets.Category, asset_cat{cat_i});
%     indx(entity.assets.Value<=0) = 0;
%     
%     DamageFunID = unique(entity.assets.DamageFunID(indx));
%     
%     for ii = 1:numel(DamageFunID)
%         fprintf('Asset DamageFunID: %d \n',DamageFunID(ii))
%         indxx = find(entity.damagefunctions.DamageFunID == DamageFunID(ii));
%         indxx = indxx(end);
%         fprintf('DamageFunID: %d, %s \n',entity.damagefunctions.DamageFunID(indxx), entity.damagefunctions.Description{indxx})
%         fprintf('max intensity %2.1f, max MDD %2.1f, \n\n', entity.damagefunctions.Intensity(indxx), entity.damagefunctions.MDD(indxx))     
%     end
% end



%%






