% create hazard plots for all hazard types, time horizons and climate
% change scenarios

%% preparations
%hazard type
% hazard_names = {'flood_depth_monsoon' 'flood_depth_cyclone' 'flood_duration_monsoon' 'flood_duration_cyclone' 'cyclone_wind'}; 
% hazard_names = {'flood_depth_monsoon' 'flood_depth_cyclone' 'flood_duration_monsoon' 'flood_duration_cyclone'}; 
hazard_names = {'flood_duration_monsoon' 'flood_duration_cyclone'}; 


%climate change scenario
cc_scenario = {'no change' 'moderate' 'extreme'}; 

%timehorizon
timehorizon = [2014 2030 2050];


%% load barisal wards
BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards_number_added.mat'];
load(BCC_wards_savename)


%% find max flood depth (monsoon and cyclone)
% hazard = barisal_hazard_entity_load('flood_depth_monsoon', cc_scenario{3}, timehorizon(2));
% event_max1 = max(full(hazard.intensity(:)));
% hazard = barisal_hazard_entity_load('flood_depth_cyclone', cc_scenario{3}, timehorizon(2));
% event_max2 = max(full(hazard.intensity(:)));
% fprintf('\n\t - Maximum flood depth from monsoon is %2.2fm.\n', event_max1)
% fprintf('\t - Maximum flood depth from cyclone is %2.2fm.\n', event_max2)
% 
% % find max flood duration (monsoon and cyclone)
% hazard = barisal_hazard_entity_load('flood_duration_monsoon', cc_scenario{3}, timehorizon(2));
% event_max1 = max(full(hazard.intensity(:)));
% hazard = barisal_hazard_entity_load('flood_duration_cyclone', cc_scenario{3}, timehorizon(2));
% event_max2 = max(full(hazard.intensity(:)));
% fprintf('\n\t - Maximum flood duration from monsoon is %2.2f days.\n', event_max1)
% fprintf('\t - Maximum flood duration from cyclone is %2.2f days.\n', event_max2)
% % event_min2 = min(nonzeros(full(hazard.intensity(:))));


%% loop over hazards
for h_i = 1:length(hazard_names)
    hazard_name = hazard_names{h_i};

    switch hazard_name
        case {'flood_depth_monsoon','flood_depth_cyclone'}
            caxis_range = [0 2.3];
        case {'flood_duration_monsoon','flood_duration_cyclone'}
            caxis_range = [0 50]; %[0 120];
    end
    % loop over time horizons
    for t_i = 1:length(timehorizon)
        
        % loop over climate change scenarios
        for cc_i = 1:length(cc_scenario)
            if ~(t_i==1 & cc_i >1) & ~(t_i>1 & cc_i== 1)
                
                % load hazard and entity and create label
                [hazard, entity, label] = barisal_hazard_entity_load(hazard_name, cc_scenario{cc_i}, timehorizon(t_i));
                
                if ~isempty(hazard) & ~isempty(entity) 
 
                    %event_sum = sum(full(hazard.intensity),2);
                    %[a,sorted_i] =sort(event_sum,'descend');
                    
                    % create hazard plots
                    for i = 1:size(hazard.intensity,1)
                        %create figure
                        fig = climada_figuresize(0.6,0.8);
                        climada_hazard_plot(hazard, i,'',caxis_range);
                        %climada_hazard_plot_hr(hazard, i);
                        %set(gca,'layer','top')
                        set(get(colorbar,'ylabel'),'string',sprintf('hazard intensity (%s)',hazard.units),'fontsize',12)
                            hazard_name_string = strrep(regexprep(hazard_name,'(\<[a-z])','${upper($1)}'),'_',' ');
                        if t_i>1
                            titlestr = sprintf('%s, (%d %s climate change), event %d', hazard_name_string,timehorizon(t_i),cc_scenario{cc_i},i);
                        else
                            titlestr = sprintf('%s (%d %s), event %d', hazard_name_string,timehorizon(t_i),cc_scenario{cc_i},i);
                        end
                        title(titlestr)
                        hold on
                        % loop over all words to plot color according to flood damage
                        BCC_ward_no  = [BCC_wards.Ward_no];
                        for ward_i = 1:length(BCC_ward_no) 
                            %indx   = find(values(ward_i)<=range_values);
                            %indx   = indx(1);
                            %indx_w = find(BCC_ward_no == ward_i);
                            plot(BCC_wards(ward_i).lon,BCC_wards(ward_i).lat,'color',[89 89 89]/255,'linewidth',1); %dark grey
                            %text(mean(BCC_wards(ward_i).lon),mean(BCC_wards(ward_i).lat),BCC_wards(ward_i).UNION_NAME,...
                            %    'Horizontalalignment','center','verticalalignment','bottom'); %grey
                        end
                        foldername = sprintf('%sresults%sdamage_plots%sHazard_%s_%d_%s_event_%d.pdf', filesep,filesep,filesep,...
                                             hazard_name,timehorizon(t_i),strrep(cc_scenario{cc_i},' ','_'),i);
                        print(gcf,'-dpdf',[climada_global.data_dir foldername])
                        close
                    end                    
                end %isempty
            end
        end %cc_i
    end %t_i
end %h_i



