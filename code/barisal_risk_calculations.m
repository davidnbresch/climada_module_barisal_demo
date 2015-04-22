%% barisal_risk_calculations

% - tc wind hazard: use barisal_tc_hazard_prob.m to create wind hazard
% - flood hazard: read asci-file from Ruud (Witteveen+Bos), load flood hazard

% flood_depth    = 1; flood_duration = 0; cyclone_wind   = 0;
% flood_duration = 1; flood_depth    = 0; cyclone_wind   = 0;
cyclone_wind   = 1; flood_depth    = 0; flood_duration = 0; 


%% hazard

if flood_depth == 1
    % hazard flood depth
    % asci_file = ;
    % hazard = climada_asci2hazard(asci_file);
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL'];
    load(hazard_set_file)
    
elseif flood_duration == 1
    % hazard duration (too be prepared!)
    % asci_file = ;
    % hazard = climada_asci2hazard(asci_file);
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_duration'];
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL'];
    load(hazard_set_file)

elseif cyclone_wind == 1 
    % hazard tc wind
    hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_prob'];
    load(hazard_set_file)

    % wind centroids
    centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
    load(centroids_file)   
end

% % load entity (asset portfolio)
% entity_file = [climada_global.data_dir filesep 'entities' filesep 'Barisal_BCC_1km_100.mat'];
% if exist(entity_file,'file')
%     load(entity_file)
% end


%% read ecorys entities (flood depth, flood duration and cyclone wind)
% barisal_entity_prepare


%% next time only load entities
if flood_depth ==1
    entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_depth.mat'];
    %entity_flood = entity;
    
elseif flood_duration ==1
    entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_duration.mat'];
    %entity_flood = entity;
    
elseif cyclone_wind == 1
    % tc wind
    entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.mat'];
    %entity_wind = entity;
end
load(entity_filename)



%% calculate damage today
if flood_depth ==1
    %entity = entity_flood;
    annotation_name = 'Flood today (depth)';
    
elseif flood_duration == 1
    %entity = entity_flood;
    annotation_name = 'Flood today (duration)';
    
elseif cyclone_wind == 1
    % tc wind
    %entity = entity_wind;
    annotation_name = 'Cyclones today';
end

force_re_encode = 1;
silent_mode     = 0;
EDS = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
climada_EDS_DFC(EDS);



%% see damage functions for different asset categories
% asset_cat = unique(entity.assets.Category);
% for cat_i = 1:length(asset_cat)
%     fprintf('-----------\n-----------\nAsset category: %s \n-----------\n',asset_cat{cat_i})
%     indx = strcmp(entity.assets.Category, asset_cat{cat_i});
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



%% waterfall graph for today, 2030, 2050
asset_cat  = unique(entity.assets.Category(entity.assets.Value>0));
entity_ori = entity;

timehorizon = [2015 2030 2050];
EDS = []; 
for t_i = 1:length(timehorizon);
    for cat_i = length(asset_cat)+1%1:length(asset_cat)+1
        
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
            %EDS_ = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
            %EDS(cat_i) = EDS_;
            EDS(t_i) = climada_EDS_calc(entity,hazard,titlestr,force_re_encode,silent_mode);
        end
    end %cat_i    
end %t_i

climada_waterfall_graph(EDS(1), EDS(2), EDS(3), 'AED')

% at the end of calculations, overwrite with original entity again
entity = entity_ori;



%% damage calculations per time horizon
asset_cat  = unique(entity.assets.Category(entity.assets.Value>0));
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
            %EDS(t_i) = climada_EDS_calc(entity,hazard,titlestr,force_re_encode,silent_mode);
        end
    end %cat_i
    
    % create figure
    figure
    climada_EDS_DFC(EDS);
    title(titlestr)
    %climada_waterfall_graph(EDS(1), EDS(2), EDS(3), 'AED')
end %t_i

% at the end of calculations, overwrite with original entity again
entity = entity_ori;



%% damage calculations per asset category
asset_cat  = unique(entity.assets.Category(entity.assets.Value>0));
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



%%






