%% prepare barisal entity, for 100 m cells basis

%hazard type
hazard_names = {'flood' 'cyclone_wind'}; 
future_years = [2014 2030 2050];

% to combine with damage functions, etc from previous ward entity file 
entity_template_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_depth.mat'];

    
%% loop over hazards
for h_i = 1:length(hazard_names)

    % set some input files
    switch hazard_names{h_i}
        case 'flood'
            hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_2014'];
            hazard_name = 'Floods';
            % entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal.xls'];
            % entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150504_values_Barisal_100m_Flooding.xls'];
            entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 040515.xls'];
            
        case 'cyclone_wind'
            hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_2014'];
            hazard_name = 'Cyclones';
            entity_filename = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Cyclones 040515.xls'];
    end
    load(hazard_set_file)


    %% read ecorys entityies
    % entity floods 2014, 2030 and 2050
    % entity cyclones 2014, 2030 and 2050
    entity       = [];
    for i = 1:length(future_years)

        entity_temp = [];
        fprintf('\n---------------\n\t - Read entity %s %d\n',hazard_name,future_years(i))

        if future_years(i) == 2014


            entity_filename_mat = strrep(entity_filename,'.xls', '.mat');
            if exist(entity_filename_mat,'file')
                delete(entity_filename_mat)
            end
            [entity,entity_save_file] = climada_entity_read(entity_filename,hazard);

            % convert local coordinates to lat lon and reencode
            fprintf('\t- Convert to lat lon and reencode\n')
            entity.assets.X = entity.assets.lon;
            entity.assets.Y = entity.assets.lat;
            [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.X, entity.assets.Y);
            entity = climada_assets_encode(entity,hazard);

            % pre-process entity in the special Barisal case (Ecorys formatting)
            entity = barisal_entity_pre_process(entity);

        else
            sheetname = sprintf('%s_%d',hazard_name,future_years(i));
            entity_temp.assets = climada_spreadsheet_read('no',entity_filename,sheetname,1);

            % convert local coordinates to lat lon and reencode
            %fprintf('\t- Convert to lat lon and reencode\n')
            %[entity_temp.assets.lon, entity_temp.assets.lat] = utm2ll_shift(entity_temp.assets.X, entity_temp.assets.Y);
            entity_temp.assets.Value       = nan(size(entity_temp.assets.X));
            entity_temp.assets.DamageFunID = nan(size(entity_temp.assets.X));
            %entity_temp = climada_assets_encode(entity_temp,hazard);

            % pre-process entity in the special Barisal case (Ecorys formatting)
            entity_temp = barisal_entity_pre_process(entity_temp);
        end

        switch future_years(i)
            case 2030
                entity.assets.Value_2030 = entity_temp.assets.Value;
            case 2050
                entity.assets.Value_2050 = entity_temp.assets.Value;
        end 
    end
    
    % we use only the assets
    assets = entity.assets;

    % to combine with damage functions, etc from previous ward entity file 
    load(entity_template_filename)
    entity.assets = assets;

    % save
    save(entity_save_file, 'entity')
    fprintf('\t -Entity %s saved\n--------------\n',hazard_name)

    
end %h_i




%% save specific entity files for flood depth, flood duration and cyclone wind speed
% hazard type
hazard_names     = {'flood_depth' 'flood_duration' 'cyclone_wind'};     
    
% loop over hazards
for h_i = 1:length(hazard_names)
     
    switch hazard_names{h_i}
        case 'flood_depth'
            entity_filename     = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 040515.mat'];
            entity_filename_new = strrep(entity_filename,'.mat', '_flood_depth.mat');
            comment             = 'Flood depth (m)';
        case 'flood_duration'
            entity_filename     = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Flooding 040515.mat'];
            entity_filename_new = strrep(entity_filename,'.mat', '_flood_duration.mat');
            comment             = 'Flood duration (days)';
        case 'cyclone_wind'
             entity_filename    = [climada_global.data_dir filesep 'entities' filesep 'Spreadsheet 100x100 Assets at risk Cyclones 040515.mat'];
            entity_filename_new = strrep(entity_filename,'.mat', '_cyclone_windspeed.mat');
            comment             = 'Cyclone wind speed (m/s)';           
    end
    load(entity_filename)
    
    % find unique damage functions and their peril_ID (flood depth, flood duration, cyclone wind speed)
    [DamFunID indx_depth indx_duration indx_windspeed] = barisal_dmgfun_filter(entity);
    
    switch hazard_names{h_i}
        case 'flood_depth'
            non_valid_DamFun    = DamFunID(~indx_depth);
        case 'flood_duration'
            non_valid_DamFun    = DamFunID(~indx_duration);
        case 'cyclone_wind'
            non_valid_DamFun    = DamFunID(~indx_windspeed);
    end
    
    if strcmp(hazard_names{h_i},'cyclone_wind')
        % check that cyclone wind speed is in m/s
        valid_indx = ~ismember(entity.damagefunctions.DamageFunID,non_valid_DamFun);
        
        % transform from kilometers per hour (kph) to m/s
        if max(entity.damagefunctions.Intensity(valid_indx)) == 400
            entity.damagefunctions.Intensity_ori         = entity.damagefunctions.Intensity;
            entity.damagefunctions.Intensity(valid_indx) = entity.damagefunctions.Intensity_ori(valid_indx)/3.6;
            entity.damagefunctions.comment               = 'Transformed intensity to m/s from km/h';
        end
    end

    % set non_hazard_type asset values to zero
    non_valid_indx                           = ismember(entity.assets.DamageFunID, non_valid_DamFun);
    entity.assets.Value(non_valid_indx)      = 0;
    entity.assets.Value_2030(non_valid_indx) = 0;
    entity.assets.Value_2050(non_valid_indx) = 0;
    entity.assets.comment                    = comment;
    % save final entity for a specific hazard type
    save(entity_filename_new, 'entity')
    fprintf('\t- Save entity %s as \n\t%s\n\n', hazard_names{h_i},entity_filename_new)
end





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




% % flood depth
% % find all assets that do not correspond to the specific index (damage function unit)
% fprintf('\t- DamageFunctions for flood depth: %d\n', numel(DamFunID(indx_depth)))
% fprintf('%d, ', DamFunID(indx_depth))
% fprintf('\n')
% entity           = entity_ori;
% non_valid_DamFun = DamFunID(~indx_depth);
% non_valid_indx   = ismember(entity.assets.DamageFunID, non_valid_DamFun);
% entity.assets.Value(non_valid_indx)      = 0;
% entity.assets.Value_2030(non_valid_indx) = 0;
% entity.assets.Value_2050(non_valid_indx) = 0;
% entity.assets.comment = 'Flood depth (m)';
% entity_filename       = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_depth.mat'];
% save(entity_filename, 'entity')
% fprintf('\t- Save entity flood depth as \n\t%s\n\n', entity_filename)
% 
% % flood duration
% % find all assets that do not correspond to the specific index (damage function unit)
% fprintf('\t- DamageFunctions for flood duration: %d\n', numel(DamFunID(indx_duration)))
% fprintf('%d, ', DamFunID(indx_duration))
% fprintf('\n')
% entity           = entity_ori;
% non_valid_DamFun = DamFunID(~indx_duration);
% non_valid_indx   = ismember(entity.assets.DamageFunID, non_valid_DamFun);
% entity.assets.Value(non_valid_indx)      = 0;
% entity.assets.Value_2030(non_valid_indx) = 0;
% entity.assets.Value_2050(non_valid_indx) = 0;
% entity.assets.comment = 'Flood duration (days)';
% entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_duration.mat'];
% save(entity_filename, 'entity')
% fprintf('\t- Save entity flood duration as \n\t%s\n\n', entity_filename)



%% hazard
% if flood == 1
%     % hazard flood
%     % asci_file = ;
%     % hazard = climada_asci2hazard(asci_file);
%     hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_2014'];
%     load(hazard_set_file)
%     hazard_name = 'Floods';
% 
% else    
%     % hazard tc wind
%     hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_2014'];
%     load(hazard_set_file)
%     hazard_name = 'Cyclones';
% 
%     % wind centroids
%     centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
%     load(centroids_file)
% end


%% --------------------
%  CYCLONE WIND
%----------------------


% %% hazard tc wind
% hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_prob'];
% load(hazard_set_file)
% 
% % wind centroids
% centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
% load(centroids_file)
  

%% read ecorys entity cyclones
% fprintf('Read entity cyclone wind\n')
% clear entity
% entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.xls'];
% entity_filename_mat = strrep(entity_filename,'.xls', '.mat');
% if exist(entity_filename_mat,'file')
%     delete(entity_filename_mat)
% end
% [entity,entity_save_file] = climada_entity_read(entity_filename,hazard);
% % convert local coordinates to lat lon
% fprintf('\t- Convert to lat lon and reencode\n')
% entity.assets.X = entity.assets.lon;
% entity.assets.Y = entity.assets.lat;
% [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.X, entity.assets.Y);
% entity = climada_assets_encode(entity,hazard);
% % save(entity_save_file, 'entity')
% % plot for first visual check
% figure
% climada_entity_plot(entity,8)
% 
% % country_name = 'Barisal';
% % check_printplot = 0;
% % printname = '';
% % keep_boundary = 0;
% % figure
% % climada_plot_entity_assets(entity,centroids,country_name,check_printplot,printname,keep_boundary);


%% organize damage functions (flood depth, flood duration, cyclone wind speed)
% %  see above
% %  find all assets that do not correspond to the specific index (damage function unit)
% fprintf('\t- DamageFunctions for wind speed: %d\n', numel(DamFunID(indx_windspeed)))
% fprintf('%d, ', DamFunID(indx_windspeed))
% fprintf('\n')
% non_valid_DamFun = DamFunID(~indx_windspeed);
% non_valid_indx   = ismember(entity.assets.DamageFunID, non_valid_DamFun);
% entity.assets.Value(non_valid_indx)      = 0;
% entity.assets.Value_2030(non_valid_indx) = 0;
% entity.assets.Value_2050(non_valid_indx) = 0;
% 
% % transform from kilometers per hour (kph) to m/s
% entity.damagefunctions.Intensity_ori = entity.damagefunctions.Intensity;
% entity.damagefunctions.Intensity     = entity.damagefunctions.Intensity_ori/3.6;
% entity.damagefunctions.comment       = 'Transformed intensity to m/s from km/h';
% 
% entity.assets.comment = 'Cyclone wind speed (m/s)';
% entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.mat'];
% save(entity_filename, 'entity')
% fprintf('\t- Save entity cyclone wind as \n\t%s\n\n', entity_filename)



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



%%






