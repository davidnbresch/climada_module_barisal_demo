%% barisal_risk_calculations

% - tc wind hazard: use barisal_tc_hazard_prob.m to create wind hazard
% - flood hazard: read asci-file from Ruud (Witteveen+Bos), load flood hazard
% prepare entity with 
% barisal_entity_prepare_100m_cells.m


%% set scenario

%hazard type
hazard_names = {'flood_depth_monsoon' 'flood_depth_cyclone' 'flood_duration_monsoon' 'flood_duration_cyclone' 'cyclone_wind'}; 

%climate change scenario
cc_scenario = {'no change' 'moderate' 'extreme'}; 

%timehorizon
timehorizon = [2014 2030 2050];


%% load barisal specifics
BCC_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_border.mat'];
load(BCC_savename)
% load BCC ward boundaries (35 polygons)
% BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards.mat'];
% BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards_Ward_no_added.mat'];
BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards_number_added.mat'];
load(BCC_wards_savename)
indx = strfind(BCC_savename,filesep);
indx2 = strfind(BCC_wards_savename,filesep);
fprintf('\t - loaded BCC specifics: %s and %s\n', BCC_savename(indx(end)+1:end), BCC_wards_savename(indx2(end)+1:end))





%% load hazard and entity
counter    = 0;
EDS_all    = '';
entity_all = '';

% loop over hazards
for h_i = 1:length(hazard_names)
    hazard_name = hazard_names{h_i};
    
    counter     = 0;
    EDS         = []; 
    silent_mode = 1;

    % loop over time horizons
    for t_i = 1:length(timehorizon)
        
        % loop over climate change scenarios
        for cc_i = 1:length(cc_scenario)
            if ~(t_i == 1 & cc_i >1)  
                
                %% load hazard and entity and create label
                [hazard, entity, label] = barisal_hazard_entity_load(hazard_name, cc_scenario{cc_i}, timehorizon(t_i));
                
                
                if ~isempty(hazard) & ~isempty(entity) 
                    
                    % loop over asset categories
                    asset_cat = unique(entity.assets.Category(entity.assets.Value>0));
                    for cat_i = length(asset_cat)+1%1:length(asset_cat)+1

                        %% the most inner loop, set counter
                        counter = counter+1;

                        % select only assets in specific category
                        if cat_i<=length(asset_cat)
                            %indx = strcmp(entity.assets.Category, asset_cat{cat_i});
                            indx = ismember(entity.assets.Category, asset_cat(1:cat_i));
                            annotation_name = asset_cat{cat_i};
                        else
                            indx = ones(size(entity.assets.Category));
                            annotation_name = 'All asset categories';
                        end
                    
                        %% calculate expected damage set (EDS)
                        entity.assets.Value(~indx) = 0;
                        force_re_encode = 0;
                        titlestr = sprintf('%s, %s, %s',label.hazard_name, label.cc_scenario, label.titlestr);
                        if isempty(EDS)
                            EDS = climada_EDS_calc(entity,hazard,titlestr,force_re_encode,silent_mode);
                            EDS.reference_year = timehorizon(t_i);
                        else
                            %EDS_ = climada_EDS_calc(entity,hazard,annotation_name,force_re_encode,silent_mode);
                            %EDS(cat_i) = EDS_;
                            EDS(counter) = climada_EDS_calc(entity,hazard,titlestr,force_re_encode,silent_mode);
                            EDS(counter).reference_year = timehorizon(t_i);
                        end
                    end %cat_i
                end %~isempty(hazard)
            end
        end %cc_i 
    end %t_i
    
    if isempty(EDS_all)
        EDS_all    = EDS;
        %entity_all = entity;
    else
        EDS_all(end+1:end+counter)    = EDS;
        %entity_all(end+1) = entity;
    end
    %fig = climada_waterfall_graph_barisal(EDS,'AED');
    %title(sprintf('Hazard type: %s',strrep(hazard_name,'_',' ')))
    %fprintf('%s: %d scenarios\n----------------\n----------------\n', upper(hazard_name), counter)
    %foldername = sprintf('%sresults%sED_waterfall_from_%s.pdf', filesep,filesep,hazard_name);
    %print(fig,'-dpdf',[climada_global.data_dir foldername])
    
end %h_i


%% write out xls
% % Percentage_Of_Value_Flag= 0;
xls_file = [climada_global.data_dir filesep 'results' filesep datestr(now,'YYYYmmDD') '_calc_ED_at_centroid_' int2str(length(EDS_all)) '_scenarios.xlsx'];
% climada_EDS_ED_at_centroid_report(EDS_all,Percentage_Of_Value_Flag,report_file)
climada_EDS_ED_at_centroid_report_xls(EDS_all,entity,xls_file)
    

 
 
%% figure for damage per ward

event_selection = [1 6 11 16 21];
t_i = 1;
% t_i = 3;
for e_i = 1:length(event_selection);
    switch e_i
        case 1
            hazard_name = 'flood depth monsoon';
        case 2
            hazard_name = 'flood depth cyclone';
        case 3
            hazard_name = 'flood duration monsoon';
        case 4
            hazard_name = 'flood duration cyclone';
        case 5
            hazard_name = 'cyclone windspeed';
    end
    
    if e_i<length(event_selection)
        fig = climada_waterfall_graph_barisal(EDS_all(event_selection(e_i):event_selection(e_i+1)-1),'AED');
    else
        fig = climada_waterfall_graph_barisal(EDS_all(event_selection(e_i):end),'AED');
    end
    title(sprintf('Hazard type: %s',strrep(hazard_name,'_',' ')))
    foldername = sprintf('%sresults%sdamage_plots%sED_waterfall_from_%s.pdf', filesep,filesep,filesep,hazard_name);
    print(fig,'-dpdf',[climada_global.data_dir foldername])
    
    %fig = climada_ED_plot_per_ward(EDS_all(event_selection(e_i)),entity,BCC_wards, timehorizon(t_i), hazard_name);
    %foldername = sprintf('%sresults%sdamage_plots%sDamage_from_%s_%d.pdf', filesep,filesep,filesep,hazard_name,timehorizon(t_i));
    %print(fig,'-dpdf',[climada_global.data_dir foldername])
    
    %fig = climada_ED_plot_per_point(EDS_all(event_selection(e_i)),entity,BCC_wards, timehorizon(t_i), hazard_name);
    %foldername = sprintf('%sresults%sdamage_plots%sDamage_per_point_from_%s_%d.pdf', filesep,filesep,filesep,hazard_name,timehorizon(t_i));
    %print(fig,'-dpdf',[climada_global.data_dir foldername])
    
    %%--ASSETS--
    %fig = climada_assets_plot_per_point(EDS_all(event_selection(e_i)+0),entity,BCC_wards, timehorizon(t_i), hazard_name);
    %foldername = sprintf('%sresults%sdamage_plots%sValue_per_point_from_%s_%d.pdf', filesep,filesep,filesep,hazard_name,timehorizon(t_i));
    %print(fig,'-dpdf',[climada_global.data_dir foldername])
    %%close
end

%%

% [hazard, entity, label] = barisal_hazard_entity_load('flood_depth_cyclone', cc_scenario{cc_i}, timehorizon(t_i));
% % climada_hazard_plot(hazard, 8)
% for i= 1:3
%     figure
%     climada_hazard_plot_hr(hazard, i);
% end
% close all


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






