%% barisal_risk_calculations

% - tc wind hazard: use barisal_tc_hazard_prob.m to create wind hazard
% - flood hazard: read asci-file from Ruud (Witteveen+Bos), load flood hazard

flood = 1;


%% hazard

if flood ~=1
    % hazard tc wind
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_prob'];
    load(hazard_set_file)

    % wind centroids
    centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
    load(centroids_file)
else
    % hazard flood
    % asci_file = ;
    % hazard = climada_asci2hazard(asci_file);
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL'];
    load(hazard_set_file)
end

% % load entity (asset portfolio)
% entity_file = [climada_global.data_dir filesep 'entities' filesep 'Barisal_BCC_1km_100.mat'];
% if exist(entity_file,'file')
%     load(entity_file)
% end

%% read ecorys entity Flooding
% entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal.xls'];
% [entity,entity_save_file] = climada_entity_read(entity_filename,hazard);
% % convert local coordinates to lat lon
% entity.assets.X = entity.assets.lon;
% entity.assets.Y = entity.assets.lat;
% [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.X, entity.assets.Y);
% entity = climada_assets_encode(entity,hazard);
% save(entity_save_file, 'entity')
% % plot for first visual check
% figure
% climada_entity_plot(entity,8)




%% read ecorys entity cyclones
% entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.xls'];
% [entity,entity_save_file] = climada_entity_read(entity_filename,hazard);
% % convert local coordinates to lat lon
% entity.assets.X = entity.assets.lon;
% entity.assets.Y = entity.assets.lat;
% [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.X, entity.assets.Y);
% entity = climada_assets_encode(entity,hazard);
% save(entity_save_file, 'entity')
% % plot for first visual check
% figure
% climada_entity_plot(entity,8)

% country_name = 'Barisal';
% check_printplot = 0;
% printname = '';
% keep_boundary = 0;
% figure
% climada_plot_entity_assets(entity,centroids,country_name,check_printplot,printname,keep_boundary);


%% organize damage functions (flood depth, flood duration, cyclone wind speed)

% find damagefunctions for flood duration (max intensity value is 14 days)
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

% flood depth
% find all assets that do not correspond to the specific index (damage function unit)
entity           = entity_ori;
non_valid_DamFun = DamFunID(~indx_depth);
non_valid_indx   = ismember(entity.assets.DamageFunID, non_valid_DamFun);
entity.assets.Value(non_valid_indx)      = 0;
entity.assets.Value_2030(non_valid_indx) = 0;
entity.assets.Value_2050(non_valid_indx) = 0;
entity.assets.comment = 'Flood depth (m)';
entity_filename       = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_depth.mat'];
save(entity_filename, 'entity')

% flood duration
% find all assets that do not correspond to the specific index (damage function unit)
entity           = entity_ori;
non_valid_DamFun = DamFunID(~indx_duration);
non_valid_indx   = ismember(entity.assets.DamageFunID, non_valid_DamFun);
entity.assets.Value(non_valid_indx)      = 0;
entity.assets.Value_2030(non_valid_indx) = 0;
entity.assets.Value_2050(non_valid_indx) = 0;
entity.assets.comment = 'Flood duration (days)';
entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_duration.mat'];
save(entity_filename, 'entity')

% cyclone wind speed
% find all assets that do not correspond to the specific index (damage function unit)
entity           = entity_ori;
non_valid_DamFun = DamFunID(~indx_duration);
non_valid_indx   = ismember(entity.assets.DamageFunID, non_valid_DamFun);
entity.assets.Value(non_valid_indx)      = 0;
entity.assets.Value_2030(non_valid_indx) = 0;
entity.assets.Value_2050(non_valid_indx) = 0;
entity.assets.comment = 'Cyclone wind speed (kph)';
entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_duration.mat'];
save(entity_filename, 'entity')





%% next time only load entities
if flood ==1
    entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal.mat'];
    load(entity_filename)
    entity_flood = entity;
else
    % tc wind
    entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.mat'];
    load(entity_filename)
    entity_wind = entity;
end


%% calculate damage today
if flood ==1
    entity = entity_flood;
    annotation_name = 'Flood today';
else
    % tc wind
    entity = entity_wind;
    annotation_name = 'Cyclones today';
end
force_re_encode = 1;
silent_mode     = 0;
EDS = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
climada_EDS_DFC(EDS);



%% see damage functions for different asset categories
asset_cat = unique(entity.assets.Category);
for cat_i = 1:length(asset_cat)
    fprintf('-----------\n-----------\nAsset category: %s \n-----------\n',asset_cat{cat_i})
    indx = strcmp(entity.assets.Category, asset_cat{cat_i});
    DamageFunID = unique(entity.assets.DamageFunID(indx));
    
    for ii = 1:numel(DamageFunID)
        fprintf('Asset DamageFunID: %d \n',DamageFunID(ii))
        indxx = find(entity.damagefunctions.DamageFunID == DamageFunID(ii));
        indxx = indxx(end);
        fprintf('DamageFunID: %d, %s \n',entity.damagefunctions.DamageFunID(indxx), entity.damagefunctions.Description{indxx})
        fprintf('max intensity %2.1f, max MDD %2.1f, \n\n', entity.damagefunctions.Intensity(indxx), entity.damagefunctions.MDD(indxx))     
    end
end


%% damage calculations per asset category
asset_cat  = unique(entity.assets.Category);
entity_ori = entity;

for cat_i = 1:length(asset_cat)+1
    
    EDS = [];
    
    % select only assets in specific category
    if cat_i<=length(asset_cat)
        indx = strcmp(entity.assets.Category, asset_cat{cat_i});
    else
        indx = ones(size(entity.assets.Category));
    end

    % risk today
    entity.assets.Value = entity_ori.assets.Value;
    entity.assets.Value(~indx) = 0;
    annotation_name = 'Risk today';
    force_re_encode = 0;
    EDS = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
    EDS(1) = EDS;
    
    % risk 2030
    entity.assets.Value = entity.assets.Value_2030;
    entity.assets.Value(~indx) = 0;
    annotation_name = 'Socio-economic 2030 (scenario 1)';
    EDS(2) = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
    
    % risk 2050
    entity.assets.Value = entity.assets.Value_2050;
    entity.assets.Value(~indx) = 0;
    annotation_name = 'Socio-economic 2050 (scenario 1)';
    EDS(3) = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
    
    % create figure
    figure
    climada_EDS_DFC(EDS);
    if cat_i<=length(asset_cat)
        title(asset_cat{cat_i})
    else
        title('All asset categories')
        % at the end of calculations, overwrite with original entity again
        entity = entity_ori;
    end
    %climada_waterfall_graph(EDS(1), EDS(2), EDS(3), 'AED')
end



%% damage calculations per time horizon
asset_cat  = unique(entity.assets.Category);
entity_ori = entity;

timehorizon = [2015 2030 2050];

for t_i = 1:length(timehorizon);
        
    EDS = [];   

    for cat_i = 1:length(asset_cat)+1
        
        switch t_i
            case 1
                % risk today
                entity.assets.Value = entity_ori.assets.Value;
                titlestr = 'Risk today';
            case 2
                % risk 2030
                entity.assets.Value = entity.assets.Value_2030;
                titlestr = 'Socio-economic 2030 (scenario 1)';
            case 3
                % risk 2050
                entity.assets.Value = entity.assets.Value_2050;
                titlestr = 'Socio-economic 2050 (scenario 1)';
        end
    
        % select only assets in specific category
        if cat_i<=length(asset_cat)
            %indx = strcmp(entity.assets.Category, asset_cat{cat_i});
            indx = ismember(entity.assets.Category, asset_cat(1:cat_i));
            annotation_name = asset_cat{cat_i};
        else
            indx = ones(size(entity.assets.Category));
            annotation_name = 'All asset categories';
        end

        entity.assets.Value(~indx) = 0;
        force_re_encode = 0;
        if isempty(EDS)
            EDS = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
        else
            EDS_ = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
            EDS(cat_i) = EDS_;
        end
    end %cat_i
    
    % create figure
    figure
    climada_EDS_DFC(EDS);
    title(titlestr)
    %climada_waterfall_graph(EDS(1), EDS(2), EDS(3), 'AED')
    
        % at the end of calculations, overwrite with original entity again
        %entity = entity_ori;
    
end %t_i




%%






