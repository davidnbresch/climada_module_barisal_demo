%% barisal_risk_calculations

% - tc wind hazard: use barisal_tc_hazard_prob.m to create wind hazard
% - flood hazard: read asci-file from Ruud (Witteveen+Bos), load flood hazard

%% set hazard type
hazard_names = {'flood_depth' 'flood_duration' 'cyclone_wind'}; 


%% load barisal specifics
BCC_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_border.mat'];
load(BCC_savename)
% load BCC ward boundaries (30 polygons)
BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards.mat'];
load(BCC_wards_savename)
indx = strfind(BCC_savename,filesep);
indx2 = strfind(BCC_wards_savename,filesep);
fprintf('\t - loaded BCC specifics: %s and %s\n', BCC_savename(indx(end)+1:end), BCC_wards_savename(indx2(end)+1:end))


%%
for h_i = 1:length(hazard_names)
    
    hazard_name = hazard_names{h_i};

    %% hazard
    switch hazard_name
        case 'flood_depth'
            % hazard = climada_asci2hazard(asci_file);
            hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL'];
            flood_depth = 1; flood_duration = 0; cyclone_wind   = 0;

        case 'flood_duration' % hazard duration (too be prepared!)
            hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_duration'];
            flood_duration = 1; flood_depth    = 0; cyclone_wind   = 0;

        case 'cyclone_wind'
            hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_prob'];
            cyclone_wind   = 1; flood_depth    = 0; flood_duration = 0; 
            % wind centroids
            %centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
            %load(centroids_file)   
    end
    load(hazard_set_file)
    fprintf('\n-----------\nSELECTED HAZARD TYPE: %s\n-----------\n', hazard_name)
    
    
    %% next time only load entities
    if flood_depth ==1
        entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_depth.mat'];

    elseif flood_duration ==1
        entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_flood_duration.mat'];

    elseif cyclone_wind == 1
        % tc wind
        entity_filename = [climada_global.data_dir filesep 'entities' filesep '20150416_values_Barisal_cyclones.mat'];
    end
    load(entity_filename)
    % create copy
    entity_ori = entity;
    indx = strfind(entity_filename,filesep);
    fprintf('\t - Loaded entity: %s\n', entity_filename(indx(end)+1:end))
    
    
    %% waterfall graph for today, 2030, 2050
    asset_cat  = unique(entity.assets.Category(entity.assets.Value>0));

    timehorizon = [2015 2030 2050];
    EDS         = []; 
    silent_mode = 1;
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
                EDS = climada_EDS_calc(entity,hazard,titlestr,force_re_encode,silent_mode);
            else
                %EDS_ = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
                %EDS(cat_i) = EDS_;
                EDS(t_i) = climada_EDS_calc(entity,hazard,titlestr,force_re_encode,silent_mode);
            end
        end %cat_i    
    end %t_i

    climada_waterfall_graph_barisal(EDS(1), EDS(2), EDS(3), 'AED');

    % at the end of calculations, overwrite with original entity again
    entity = entity_ori;
    
    
    
    %% write annual damage results to excel
    buffer_hazard = 0;
    if flood_depth == 1
        indx = strfind(entity_filename,'_');
        xls_filename = sprintf('%s.xls',entity_filename(1:indx(end-1)-1));
        sheet_name = 'assets';
        
    elseif flood_duration == 1   
        indx = strfind(entity_filename,'_');
        xls_filename = sprintf('%s.xls',entity_filename(1:indx(end-1)-1));
        sheet_name = 'assets';
        buffer_hazard = 4;
        
    elseif cyclone_wind == 1
        % tc wind
        indx = strfind(entity_filename,'_');
        xls_filename = sprintf('%s.xls',entity_filename(1:indx(end)-1));
        sheet_name = 'assets cyclones';
        buffer_hazard = 0;
    end
    filename_add = sprintf('_annual_damage_calculated_%s.xls',datestr(now,'YYYYmmDD'));
    new_xls_filename = strrep(xls_filename,'.xls',filename_add);
    new_xls_filename = strrep(new_xls_filename,'entities','results');
    
    if ~exist(new_xls_filename,'file')
        % copy and rename xls-file
        copyfile(xls_filename,new_xls_filename);
    end
    
    % add annual damage column to existing xls-file 
    existing_row_no = numel(fieldnames(entity.assets));
    buffer_rows     = 3;
    row_no          = existing_row_no+buffer_rows+buffer_hazard;
    alphabet_no     = 26;
    if row_no>=alphabet_no
        multiple = row_no/alphabet_no;
        row_no = rem(row_no, alphabet_no);
        multiple = char('A'+multiple-1);
    else
        multiple = '';
    end
    xlRange         = sprintf('%s%s1',multiple,char('A'+row_no));
    fprintf('\t - Results into xls: Write annual damages from %s to sheet "%s" (row %s)\n', strrep(hazard_name,'_',' '), sheet_name, xlRange(1:end-1))
    matr     = cell(length(EDS(1).ED_at_centroid)+3,numel(EDS));
    for e_i = 1:numel(EDS)
        matr{1,e_i} = sprintf('Annual damage from %s, %d, %s', strrep(hazard_name,'_',' '), timehorizon(e_i), EDS(e_i).annotation_name);
        matr(2:end-2,e_i) = num2cell(EDS(e_i).ED_at_centroid);
        matr{end-1,e_i}   = sprintf('Total damage');
        matr{end,e_i}     = sum(EDS(e_i).ED_at_centroid);
    end
    xlswrite(new_xls_filename, matr,sheet_name,xlRange)
    
    %% figure for damage per ward
    t_i = 1;
    fig = climada_ED_plot_per_ward(EDS(t_i),entity,BCC_wards, timehorizon(t_i), hazard_name);
    foldername = sprintf('%sresults%sDamage_from_%s_%d.pdf', filesep,filesep,hazard_name,timehorizon(t_i));
    print(fig,'-dpdf',[climada_global.data_dir foldername])
    %close
    
    
end %h_i


%% check hazard intensity and max hazard event
% [max_int, max_event] = max(full(sum(hazard.intensity,2)));
% % % [max_int, max_event] = max(full(max(hazard.intensity,[],2)));
% % % climada_hazard_footprint_plot(hazard, max_event, '');
% % hazard_sum = hazard;
% % hazard_sum.intensity(1,:) = sum(hazard_sum.intensity);
% 
% % for e_i = 1:1:29
%     fig = climada_figuresize(0.75,0.75);
%     %%climada_hazard_plot(hazard_sum,1)
%     %e_i = 27;
%     %e_i = 8;
%     e_i = max_event;
%     climada_hazard_plot(hazard,e_i);
%     
%     % plot Wards
%     hold on
%     for w_i=1:length(BCC_wards)
%         h = plot(BCC_wards(w_i).lon,BCC_wards(w_i).lat,'color',[244 164 96 ]/255);%sandybrown;
%         %text(mean(BCC_wards(w_i).lon),mean(BCC_wards(w_i).lat),sprintf('Ward %d',BCC_wards(w_i).Ward_no));%sandybrown;
%         g = text(entity.assets.lon(w_i), entity.assets.lat(w_i),sprintf('Ward %d',w_i),'Horizontalalignment','center','verticalalignment','bottom');
%     end
%     indx_a = entity.assets.Value>0;
%     g = plot(entity.assets.lon(indx_a),entity.assets.lat(indx_a),'kx','markersize',8,'LineWidth',1.5);
%     legend([h g], 'Ward boundaries', 'Lat/lon coordinates for assets')
%     xlabel('Longitude'); ylabel('Latitude')
%     foldername = sprintf('%sresults%sHazard_%s_event_%d.pdf',filesep, filesep,hazard_name,e_i);
%     print(fig,'-dpdf',[climada_global.data_dir foldername])
%     close
% % end

% %figure
% %[h h_points] = plotclr(hazard.lon,hazard.lat,hazard.intensity(e_i,:), 'o', 5, 1, 0, 0.5, '', 0, 0);
% 
% % ward11_no = 11;
% % entity.assets.Value(ward11_no);
% % plot(entity.assets.lon(ward11_no), entity.assets.lat(ward11_no),'dr')
% % c = entity.assets.centroid_index(ward11_no);
% % plot(hazard.lon(c), hazard.lat(c),'pb')
% % int = full(hazard.intensity(:,c));
% % plot(hazard.lon, hazard.lat, 'x')
% % int = full(hazard.intensity(:,c+128));
% % int = full(hazard.intensity(:,c+0));
% % 
% % ward6_no = 6;
% % plot(entity.assets.lon(ward6_no), entity.assets.lat(ward6_no),'dr')
% % c6 = entity.assets.centroid_index(ward6_no);
% % plot(hazard.lon(c6), hazard.lat(c6),'pb')
% % int = full(hazard.intensity(:,c6));
% % end
% 
% % % load entity (asset portfolio)
% % entity_file = [climada_global.data_dir filesep 'entities' filesep 'Barisal_BCC_1km_100.mat'];
% % if exist(entity_file,'file')
% %     load(entity_file)
% % end


%% read ecorys entities (flood depth, flood duration and cyclone wind)
% barisal_entity_prepare







%% calculate damage today
% if flood_depth ==1
%     annotation_name = 'Flood today (depth)';
%     
% elseif flood_duration == 1
%     annotation_name = 'Flood today (duration)';
%     
% elseif cyclone_wind == 1
%     % tc wind
%     annotation_name = 'Cyclones today';
% end
% 
% force_re_encode = 1;
% silent_mode     = 0;
% EDS = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
% climada_EDS_DFC(EDS);



%% see damage functions for different asset categories
% % asset_cat = unique(entity.assets.Category);
% asset_cat  = unique(entity.assets.Category(entity.assets.Value>0));
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


%% figure to check if max hazard (flood duration) touches upon asset values
% figure
% [max_int, max_event] = max(full(sum(hazard.intensity,2)));
% climada_hazard_footprint_plot(hazard, max_event, '');
% asset_cat  = unique(entity.assets.Category(entity.assets.Value>0));
% color_ = jet(length(asset_cat));
% marker_ = {'+','o','d','*','s','p','x'};
% for cat_i = 1:length(asset_cat)
%     indx = strcmp(entity.assets.Category,asset_cat(cat_i));
%     hold on
%     if cat_i<=length(marker_)
%         h(cat_i) = plot(entity.assets.lon(indx), entity.assets.lat(indx),marker_{cat_i},'color',color_(cat_i,:));
%     else
%         h(cat_i) = plot(entity.assets.lon(indx), entity.assets.lat(indx),'p','color',color_(cat_i,:));
%     end
% end
% legend(h,asset_cat)



%% figure BCC wards and numbers
% climada_figuresize(0.7,0.5)
% for w_i=1:length(BCC_wards)
%     hold on
%     h(3)= plot(BCC_wards(w_i).lon,BCC_wards(w_i).lat,'color',[244 164 96 ]/255);%sandybrown;
%     %text(mean(BCC_wards(w_i).lon), mean(BCC_wards(w_i).lat), int2str(BCC_wards(w_i).WARDS_F_ID))
%     text(mean(BCC_wards(w_i).lon), mean(BCC_wards(w_i).lat), int2str(BCC_wards(w_i).Ward_no))
% end
% plot(entity.assets.lon(1:30),entity.assets.lat(1:30), 'x')
% for a_i = 1:30
%     text(entity.assets.lon(a_i),entity.assets.lat(a_i), int2str(entity.assets.Ward(a_i)))
% end


%% assign ward numbers to ward shapefile
% ward_points = [entity.assets.lon(1:30) entity.assets.lat(1:30)];
% for w_i=1:length(BCC_wards)
%     polygon_nodes = [BCC_wards(w_i).lon' BCC_wards(w_i).lat'];
%     [cn,on] = inpoly(ward_points,polygon_nodes);
%     BCC_wards(w_i).Ward_no = find(cn);
% end
% % hand corrections
% BCC_wards(2).Ward_no = 6;
% BCC_wards(3).Ward_no = 5;
% BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards.mat'];
% save(BCC_wards_savename,'BCC_wards')







%% damage calculation per ward (per category, for today, 2030 and 2050)
% silent_mode= 1;
% asset_cat  = unique(entity.assets.Category(entity.assets.Value>0));
% ward_no    = unique(entity.assets.Ward);
% 
% timehorizon = [2014 2030 2050];
% EDS         = [];
% for t_i = 1:length(timehorizon);
%      
%     counter     = 0;
% 
%     % select ward
%     for ward_i = 1:length(ward_no) 
%         
%         % select only one ward
%         if ward_i<=length(ward_no)
%             indx_w = ismember(entity.assets.Ward, ward_no(ward_i));
%             %indx_w = ismember(entity.assets.Ward, ward_no(1:ward_i));
%             annotation_name_ward = sprintf('Ward %d',ward_no(ward_i));
%         else
%             indx_w = ones(size(entity.assets.Ward));
%             annotation_name_ward = 'All Wards';
%         end
%         
%         % select asset category
%         for cat_i = 1:length(asset_cat)+1%length(asset_cat)+1 
% 
%             %select time horizone
%             switch t_i 
%                 case 1
%                     % risk today
%                     entity.assets.Value = entity_ori.assets.Value;
%                     titlestr = 'Risk today';
%                 case 2
%                     % risk 2030
%                     entity.assets.Value = entity.assets.Value_2030;
%                     titlestr = 'Socio-economic 2030 (scenario 1)';
%                 case 3
%                     % risk 2050
%                     entity.assets.Value = entity.assets.Value_2050;
%                     titlestr = 'Socio-economic 2050 (scenario 1)';
%             end
% 
%             % select only assets in specific category
%             if cat_i<=length(asset_cat)
%                 %indx = strcmp(entity.assets.Category, asset_cat{cat_i});
%                 indx = ismember(entity.assets.Category, asset_cat(1:cat_i));
%                 annotation_name = asset_cat{cat_i};
%             else
%                 indx = ones(size(entity.assets.Category));
%                 annotation_name = 'All asset categories';
%             end
%             
%             %------
%             % combine ward and category index
%             %------
%             indx      = logical(indx);
%             indx_w    = logical(indx_w);
%             indx_comb = logical(indx.*indx_w);
%             counter   = counter+1;
%             annotation_name_comb = sprintf('%d, %s, %s', timehorizon(t_i), annotation_name_ward, annotation_name);
%             
%             % a quick check
%             %entity.assets.Ward(indx_comb)
%             %entity.assets.Category(indx_comb)
%             %entity.assets.Value(indx_comb)
% 
%             entity.assets.Value(~indx_comb) = 0;
%             force_re_encode = 0;
%             if isempty(EDS)
%                 EDS = climada_EDS_calc(entity,hazard,annotation_name_comb,force_re_encode,silent_mode);
%             else
%                 %EDS_ = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
%                 %EDS(cat_i) = EDS_;
%                 EDS(counter) = climada_EDS_calc(entity,hazard,annotation_name_comb,force_re_encode,silent_mode);
%             end
%             
%         end %cat_i  
%     end %ward_i
% end %t_i
% % climada_waterfall_graph_barisal(EDS(1), EDS(2), EDS(3), 'AED')
% 
% % at the end of calculations, overwrite with original entity again
% entity = entity_ori;
% 
% return

%% figures for specific asset categories and time horizons
% EDS_annotation_names = {EDS.annotation_name};
% for t_i = 1:length(timehorizon);
%     for cat_i = 1:length(asset_cat)+1
%         
%         % select only assets in specific category
%         if cat_i<=length(asset_cat)           
%             annotation_name = asset_cat{cat_i};
%         else
%             annotation_name = 'All asset categories';
%         end
%         indx_2 = strfind(EDS_annotation_names, annotation_name);
%         indx_2 = find(~cellfun(@isempty,indx_2)); 
%         
%         indx_3 = strfind(EDS_annotation_names, sprintf('%d',timehorizon(t_i)));
%         
%         %-------
%         %annual expected damage
%         %-------
%         fig = climada_figuresize(0.75,0.75);
%         no_colors    = 10;
%         cbar         = jet(no_colors);
%         ED           = [EDS(indx_2).ED];
%         pos_indx     = ED>0;
%         min_value    = min(ED(pos_indx));
%         max_value    = max(ED);
%         range_values = linspace(min_value,max_value,no_colors);
%         BCC_ward_no  = [BCC_wards.Ward_no];
%         for ward_i = 1:length(ward_no) 
%             indx = find(EDS(indx_2(ward_i)).ED<=range_values);
%             indx = indx(1);
%             %h(3)= plot(BCC_wards(w_i).lon,BCC_wards(w_i).lat,'color',[244 164 96 ]/255);%sandybrown
%             indx_w = find(BCC_ward_no == ward_i);
%             h(3)   = fill(BCC_wards(indx_w).lon,BCC_wards(indx_w).lat,cbar(indx,:));
%             hold on
%             indx_a = find(EDS(indx_2(ward_i)).assets.Value);
%             if ~isempty(indx_a); indx_a = indx_a(1);end
%             indx_text = strfind(EDS(indx_2(ward_i)).annotation_name,',');
%             if ~isempty(indx_a)
%                 text(EDS(indx_2(ward_i)).assets.lon(indx_a), EDS(indx_2(ward_i)).assets.lat(indx_a), EDS(indx_2(ward_i)).annotation_name(1:indx_text-1),...
%                     'Horizontalalignment','center','verticalalignment','bottom','color','k') %brighten(cbar(indx,:),-0.8)
%                 g = plot(mean(EDS(indx_2(ward_i)).assets.lon(indx_a)), mean(EDS(indx_2(ward_i)).assets.lat(indx_a)),'kx','markersize',8,'LineWidth',1.5);
%             end
%         end
%         colormap(cbar)
%         t = colorbar;
%         %cbar_label = sprintf('Intensity %s (%s)', hazard.peril_ID, hazard.units);
%         set(get(t,'ylabel'),'String', ('1000 BDT'),'fontsize',12);
%         caxis([min_value max_value])
%         %axislim = [min(EDS(1).assets.lon) max(EDS(1).assets.lon)*1 min(EDS(1).assets.lat) max(EDS(1).assets.lat)*1];
%         %axislim = [min(hazard.lon) max(hazard.lon)*1 min(hazard.lat) max(hazard.lat)*1];
%         %axislim = [90.25 90.45 22.6 22.8]; %barisal close up BCC 
%         axislim = [90.297 90.3957 22.64 22.752]; %barisal close up BCC 
%         axis(axislim)
%         axis equal
%         %titlestr = sprintf('%d, Annual damage, %s - %s', timehorizon(t_i), EDS(1).annotation_name, EDS(end).annotation_name);
%         titlestr = sprintf('%d: Annual damage from %s', timehorizon(t_i), strrep(hazard_name,'_',' '));
%         title({titlestr; annotation_name})
%         legend(g,'Lat/lon coordinates for assets')
%         foldername = sprintf('%sresults%sDamage_from_%s_%d.pdf', filesep,filesep,hazard_name,timehorizon(t_i));
%         print(fig,'-dpdf',[climada_global.data_dir foldername])
%         %close
% 
%         %-------
%         %values
%         %-------
% %         fig = climada_figuresize(0.75,0.75);
% %         no_colors    = 10;
% %         cbar         = jet(no_colors);
% %         Value        = [EDS.Value];
% %         pos_indx     = Value>0;
% %         min_value    = min(Value(pos_indx));
% %         max_value    = max(Value);
% %         range_values = linspace(min_value,max_value,no_colors);
% %         BCC_ward_no  = [BCC_wards.Ward_no];
% %         for ward_i = 1:length(ward_no) 
% %             indx   = find(EDS(ward_i).Value<=range_values);
% %             indx   = indx(1);
% %             %h(3)  = plot(BCC_wards(w_i).lon,BCC_wards(w_i).lat,'color',[244 164 96 ]/255);%sandybrown
% %             indx_w = find(BCC_ward_no == ward_i);
% %             h(3)   = fill(BCC_wards(indx_w).lon,BCC_wards(indx_w).lat,cbar(indx,:));
% %             indx_a = find(EDS(ward_i).assets.Value);
% %             indx_a = indx_a(1);
% %             indx_text = strfind(EDS(ward_i).annotation_name,',');
% %             text(EDS(ward_i).assets.lon(indx_a), EDS(ward_i).assets.lat(indx_a), EDS(ward_i).annotation_name(1:indx_text-1),...
% %                 'Horizontalalignment','center','verticalalignment','bottom','color','k')%brighten(cbar(indx,:),-0.8
% %             hold on
% %             g = plot(mean(EDS(ward_i).assets.lon(indx_a)), mean(EDS(ward_i).assets.lat(indx_a)),'kx','markersize',8,'LineWidth',1.5);
% %         end
% %         colormap(cbar)
% %         t = colorbar;
% %         %cbar_label = sprintf('Intensity %s (%s)', hazard.peril_ID, hazard.units);
% %         set(get(t,'ylabel'),'String', ('1000 BDT'),'fontsize',12);
% %         caxis([min_value max_value])
% %         axis(axislim)
% %         axis equal
% %         %titlestr = sprintf('%d, Values, %s - %s', timehorizon(t_i), EDS(1).annotation_name, EDS(end).annotation_name);
% %         titlestr = sprintf('%d: Values', timehorizon(t_i));
% %         title({titlestr; annotation_name})
% %         legend(g,'Lat/lon coordinates for assets')
% %         foldername = sprintf('%sresults%sValues_for_%s_%d.pdf', filesep,filesep,hazard_name,timehorizon(t_i));
% %         print(fig,'-dpdf',[climada_global.data_dir foldername])
% %         close
%     end %cat_i
% end %t_i
    




%% damage calculations per time horizon
% asset_cat  = unique(entity.assets.Category(entity.assets.Value>0));
% entity_ori = entity;
% 
% timehorizon = [2015 2030 2050];
% 
% for t_i = 2%:length(timehorizon);
%         
%     EDS = [];   
% 
%     for cat_i = 1:length(asset_cat)+1
%         
%         switch t_i
%             case 1
%                 % risk today
%                 entity.assets.Value = entity_ori.assets.Value;
%                 titlestr = 'Risk today';
%             case 2
%                 % risk 2030
%                 entity.assets.Value = entity.assets.Value_2030;
%                 titlestr = 'Socio-economic 2030 (scenario 1)';
%             case 3
%                 % risk 2050
%                 entity.assets.Value = entity.assets.Value_2050;
%                 titlestr = 'Socio-economic 2050 (scenario 1)';
%         end
%     
%         % select only assets in specific category
%         if cat_i<=length(asset_cat)
%             %indx = strcmp(entity.assets.Category, asset_cat{cat_i});
%             indx = ismember(entity.assets.Category, asset_cat(1:cat_i));
%             annotation_name = asset_cat{cat_i};
%         else
%             indx = ones(size(entity.assets.Category));
%             annotation_name = 'All asset categories';
%         end
% 
%         entity.assets.Value(~indx) = 0;
%         force_re_encode = 0;
%         if isempty(EDS)
%             EDS = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
%         else
%             EDS_ = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
%             EDS(cat_i) = EDS_;
%             %EDS(t_i) = climada_EDS_calc(entity,hazard,titlestr,force_re_encode,silent_mode);
%         end
%     end %cat_i
%     
%     % create figure
%     figure
%     climada_EDS_DFC(EDS);
%     title(titlestr)
%     %climada_waterfall_graph_barisal(EDS(1), EDS(2), EDS(3), 'AED')
% end %t_i
% 
% % at the end of calculations, overwrite with original entity again
% entity = entity_ori;



%% damage calculations per asset category
% asset_cat  = unique(entity.assets.Category(entity.assets.Value>0));
% entity_ori = entity;
% 
% for cat_i = 1:length(asset_cat)+1
%     
%     EDS = [];
%     
%     % select only assets in specific category
%     if cat_i<=length(asset_cat)
%         indx = strcmp(entity.assets.Category, asset_cat{cat_i});
%     else
%         indx = ones(size(entity.assets.Category));
%     end
% 
%     % risk today
%     entity.assets.Value = entity_ori.assets.Value;
%     entity.assets.Value(~indx) = 0;
%     annotation_name = 'Risk today';
%     force_re_encode = 0;
%     EDS = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
%     EDS(1) = EDS;
%     
%     % risk 2030
%     entity.assets.Value = entity.assets.Value_2030;
%     entity.assets.Value(~indx) = 0;
%     annotation_name = 'Socio-economic 2030 (scenario 1)';
%     EDS(2) = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
%     
%     % risk 2050
%     entity.assets.Value = entity.assets.Value_2050;
%     entity.assets.Value(~indx) = 0;
%     annotation_name = 'Socio-economic 2050 (scenario 1)';
%     EDS(3) = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
%     
%     % create figure
%     figure
%     climada_EDS_DFC(EDS);
%     if cat_i<=length(asset_cat)
%         title(asset_cat{cat_i})
%     else
%         title('All asset categories')
%         % at the end of calculations, overwrite with original entity again
%         entity = entity_ori;
%     end
%     %climada_waterfall_graph_barisal(EDS(1), EDS(2), EDS(3), 'AED')
% end



%%






