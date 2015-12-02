%% Barisal Risk Calculations

% Lea Mueller, muellele@gmail.com, 20150906, baseline run for scenarios 1, 2, 3
% Lea Mueller, muellele@gmail.com, 20151117, project package run for scenario 1 based on a 
%           - new assets (location, values and damage function change)
%           - reduced hazard (FL depth monsoon, FL duration mosoon, FL depth cyclone, FL duration cyclone)
%           - land raising by 0.3 m(hazard_intensity_impact_b = -0.3)

clc
climada_global.waitbar = 0;
climada_global.EDS_at_centroid = 0;

%% Directories
barisal_data_dir = climada_global.data_dir;
% barisal_data_dir= [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
entities_dir    = [barisal_data_dir filesep 'entities'];
hazards_dir     = [barisal_data_dir filesep 'hazards'];
results_dir     = [barisal_data_dir filesep 'results' filesep '20151116_measure_project_package'];
% results_dir     = [barisal_data_dir filesep 'results' filesep '20150909_new_runs'];
% results_dir     = [barisal_data_dir filesep 'results'];

%% load barisal specifics
BCC_border_file = [entities_dir filesep 'BCC_border.mat'];
load(BCC_border_file)

% load BCC ward boundaries (35 polygons)
BCC_wards_file = [entities_dir filesep 'BCC_wards_number_added.mat'];
load(BCC_wards_file)
BCC_wards_ll     = BCC_wards;
[BCC_wards_ll.x] = BCC_wards_ll.X;      BCC_wards_ll = rmfield(BCC_wards_ll,'X');
[BCC_wards_ll.y] = BCC_wards_ll.Y;      BCC_wards_ll = rmfield(BCC_wards_ll,'Y');
[BCC_wards_ll.X] = BCC_wards_ll.lon;    BCC_wards_ll = rmfield(BCC_wards_ll,'lon');
[BCC_wards_ll.Y] = BCC_wards_ll.lat;    BCC_wards_ll = rmfield(BCC_wards_ll,'lat');

climada_admin_name('Bangladesh','Barisal S.',4,1);

clear BCC_border_file BCC_wards_file


%% hazard file name

% This section stores the filenames of all the relevant hazards in a cell
% array for easy retrieval later (using barisal_get_hazard)

% Step 1: call barisal_hazard_read to read Witteveen & Bos hazard asci
% files and save resulting hazard structs as .mat files in the climada_GIT
% barisal_demo module data directory.
force_hazard_asci_read = 0;
if force_hazard_asci_read
    barisal_hazard_read
end

% Step 2: construct a cell array with all the hazard file names
hazard_file_tmp = 'Barisal_BCC_hazard_PIDSPEC_CCSCEN.mat';

% loop over climate change scenarios
CCSCEN  = {'2014' 'cc_2030_moderate' 'cc_2030_extreme' 'cc_2050_moderate' 'cc_2050_extreme'};
% loop over hazard perils
PID = {'TC' 'FL_depth_cyclone' 'FL_depth_monsoon' 'FL_duration_cyclone' 'FL_duration_monsoon'};
SPEC= {'' }; %'_rain_only'};

file_i = 0; hazard_files = {};
for cc = CCSCEN
    for pid = PID
        for spec = SPEC
            hazard_file = hazard_file_tmp;
            hazard_file = strrep(hazard_file,'PID',char(pid));
            hazard_file = strrep(hazard_file,'SPEC',char(spec));
            hazard_file = strrep(hazard_file,'CCSCEN',char(cc));
            if exist([hazards_dir filesep hazard_file],'file')
                file_i = file_i+1;
                hazard_files{file_i} = [hazards_dir filesep hazard_file];
                load(hazard_files{file_i})
                % ensure peril ID is correct
                if ~isempty(strfind(hazard_file,'duration'))
                    hazard.peril_ID = 'FL_';
                end
                % consistency in hazard.comment
                if strcmp(hazard.peril_ID,'TC')
                    % hazard.comment = strrep(strrep(hazard_file,'_',' '),'.mat','');
                    hazard.comment = ['TC wind speed ' strrep(char(cc),'_',' ')];
                    fprintf('%s\n',hazard.comment)
                else
                    hazard.comment = strrep([char(pid) ' ' char(cc) ', modelled by W+B'],'_',' ');
                    fprintf('%s\n',hazard.comment)
                end
                
                % datenum for monsoon hazard events
                if ~isempty(strfind(pid,'monsoon'))
                    hazard.yyyy     = [1:length(hazard.yyyy)]+1982;
                    hazard.datenum  = datenum(hazard.yyyy,hazard.mm,hazard.dd);
                end
                
                % consistency in hazard.filename
                hazard.filename = hazard_files{file_i};
                % save corrected hazard
                save(hazard.filename,'hazard')
            end
        end
    end
end
clear cc CCSCEN pid PID spec SPEC file_i hazard_file_tmp hazard_file force_hazard_asci_read

return 


%% entity

% This section constructs a cell array with entity file names for easy
% retrieval later with barisal_get_entity

% Entity excel file from Ecorys (unstructured)
% entity_file_xls = [entities_dir filesep 'BCC_entity_260615_se_1.xls'];
% entity_file_xls = [entities_dir filesep 'Asset Entity CLIMADA - Effect Project Package 07072015.xls'];
% run 20151116 measure project package
entity_file_xls{1} = [entities_dir filesep '20151116_measure_project_package' filesep 'Assets_at_risk_100x100_16112015 Final FS Project package - scenario 1.xlsx'];
% run 20150909 baseline
% entity_file_xls{1} = [entities_dir filesep '20150909_new_runs' filesep 'Assets_at_risk_100x100_09092015 Baseline scenario 1.xlsx'];
% entity_file_xls{2} = [entities_dir filesep '20150909_new_runs' filesep 'Assets_at_risk_100x100_09092015 Baseline scenario 2.xlsx'];
% entity_file_xls{3} = [entities_dir filesep '20150909_new_runs' filesep 'Assets_at_risk_100x100_09092015 Baseline scenario 3.xlsx'];

% % Different entities for the adaptation measures
% % entity_file_xls = [entities_dir filesep 'Measure_spatial_planning_290615_se_1.xls'];
% % entity_file_xls = [entities_dir filesep 'Measure_early_warning_system_290615_se_1.xls'];
% % entity_file_xls = [entities_dir filesep 'Measure_Flood_proof_road_infrastructure_290615_se_1.xls'];
% entity_file_xls = [entities_dir filesep 'Measure_package_070715_se_1.xls'];

% Damage function file from Ecorys
damfun_file_xls = [entities_dir filesep 'BCC_dmg_functions_260615.xls'];

% Entity template file from global data dir
entity_temp_xls = [climada_global.data_dir filesep 'entities' filesep 'entity_template.xls'];

% Sheet names in Ecorys entity xls file
sheets          = {'Floods_2014' 'Floods_2030' 'Floods_2050' 'Cyclones_2014' 'Cyclones_2030' 'Cyclones_2050'};
% sheets          = {'2014' '2030' '2050'}; % for population entity

% Whether to re-read entity/damagefunctions from Ecorys xls (1) or load matfiles (0)
force_assets_re_read   = 0;
force_damfun_re_read   = 0;
add_measures = 1;
counter = 0;
% counter = 6;

for scenario_eco_i = 1:numel(entity_file_xls)
    for s_i = 1:length(sheets)
        clear entity

        % mat file name - separate mat file for each sheet in Ecorys entity xls
        entity_file_mat     = strrep(entity_file_xls{scenario_eco_i},'.xlsx',['_' lower(sheets{s_i}) '.mat']);
        counter = counter+1;
        
        if exist(entity_file_mat,'file') && ~force_assets_re_read
            % Load mat file
            [~,fN] = fileparts(entity_file_mat);
            fprintf('entity %s already exists, skipping\n',fN)

            entity_files{counter} = entity_file_mat;

    %         load(entity_file_mat)
    %         % append filename for consistency
    %         entity.assets.filename = entity_file_mat;
    % 
    %         % get reference year, comment from sheet name
    %         [~,yr_]                         = strtok(sheets{s_i},'_');
    %         entity.assets.reference_year    = str2num(yr_(2:end));
    % 
    %         % add income information to entity.assets for residential categories only
    %         entity = barisal_entity_pre_process_income(entity);
    %         
    %         % find top-most vulnerable buildings -- load specific EDS before
    %         switch entity.assets.reference_year
    %             case 2030,  criterion_A = 0.25;     criterion_B = 0.10;
    %             case 2050,  criterion_A = 0.50;     criterion_B = 0.15;
    %             otherwise, continue;
    %         end
    %         load([entities_dir filesep 'ED_at_centroid_baseline.mat'])
    %         entity    = barisal_ED_find_most_vulnerable(entity, EDS_FL, criterion_A,criterion_B);
    % 
    %         % save corrected entity
    %         save(entity.assets.filename,'entity')

    % ==================
    % This section for looking at slum houses only. Be careful when
    % uncommenting!
    % ==================
    %         entity.assets.Value(~(strcmp(entity.assets.Category,'Residential_buildings_Juphri_ASSETS')...
    %             | strcmp(entity.assets.Category,'Residential_buildings_Katcha_ASSETS')))=0;
    %
    %         figure; hold on
    %         climada_entity_plot(entity,5,'BDT',0);
    %         box on
    %         shape_plotter(BCC_wards_ll,'','linewidth',1,'color',[81 81 81]/255);
    %         axis([min(entity.assets.lon)-0.01 max(entity.assets.lon)+0.01 min(entity.assets.lat)-0.01 max(entity.assets.lat)+0.01])
    %
    %         entity_file_mat = strrep(entity_file_mat,'.mat','_slums.mat');
    %         [~,fN] = fileparts(entity_file_mat);
    %         print(gcf,'-dpng',[results_dir filesep fN '.png'])
    %         clear fN
    %         close
    %====================

        else
            force_damfun_re_read = 1; % since file does not yet exist
            % read in entity from Ecorys xls
            % assets
            entity.assets = climada_xlsread('no',entity_file_xls{scenario_eco_i},sheets{s_i},1);
            % restructure into climada entity format (not necessary for
            % population entity)
            entity        = barisal_entity_pre_process(entity);
            % remove NaN entries in assets uniformly from all fields
            nan_ndx = isnan(entity.assets.Ward_Nr);
            flds = fieldnames(entity.assets);
            for fld_i = 1:length(flds)
                if (iscell(entity.assets.(flds{fld_i})) || isnumeric(entity.assets.(flds{fld_i}))) ...
                        && length(entity.assets.(flds{fld_i}))==length(nan_ndx)
                    entity.assets.(flds{fld_i})(nan_ndx) = [];
                end
            end
            entity.assets.Value(isnan(entity.assets.Value)) = 0;

            % set negative asset values to 0
            fprintf('setting %d negative asset values to zero\n',sum(entity.assets.Value<0))
            entity.assets.Value(entity.assets.Value<0) = 0;

            % set deductible & cover
            fprintf('setting deductible and cover\n')
            entity.assets.Deductible    = entity.assets.Value .* 0;
            entity.assets.Cover         = entity.assets.Value;

            % add income information to entity.assets for residential categories only
            %entity = barisal_entity_pre_process_income(entity);

            % find top-most vulnerable buildings -- load specific EDS before
            %criterion = 0.050;
            %entity    = barisal_ED_find_most_vulnerable(entity, EDS, criterion);

            % coord transformation from UTM to lat lon
            [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.lon, entity.assets.lat);

            % get reference year, comment from sheet name
            [region,yr_]                    = strtok(sheets{s_i},'_');
            entity.assets.reference_year    = str2num(yr_(2:end));
            entity.assets.comment           = strrep(sheets{s_i},'_',' ');
            entity.assets.filename          = entity_file_mat;
            entity.assets.region            = sprintf('Barisal assets, exposed to %s, Scenario %d',lower(region),scenario_eco_i);

            fprintf('saving entity as %s\n',entity.assets.filename)
            save(entity.assets.filename,'entity')
            entity_files{counter} = entity.assets.filename; % for filename consistency
        end

        if force_damfun_re_read
            fprintf('adding damage functions from %s\n',damfun_file_xls)
            load(entity_file_mat)
            
            %add region
            [region,yr_]                    = strtok(sheets{s_i},'_');
            entity.assets.region            = sprintf('Barisal assets exposed to %s',lower(region));
            
            % damagefunctions from Ecorys xls
            entity.damagefunctions          = climada_xlsread(0,damfun_file_xls,'formatted',1);

            % find TC damfuns
            tc_ndx                          = strcmp(entity.damagefunctions.peril_ID,'TC');
            % widnspeed conversion from km/h to m/s
            entity.damagefunctions.Intensity(tc_ndx) = entity.damagefunctions.Intensity(tc_ndx)./3.6;
            entity.damagefunctions.units    (tc_ndx) = {'m/s'};

            % need MDD and PAA for EDS_calc, not given in Ecorys entity, so
            % take sqrt of MDR instead.
            entity.damagefunctions.MDD = sqrt(entity.damagefunctions.MDR);
            entity.damagefunctions.PAA = sqrt(entity.damagefunctions.MDR);

            % discount
            entity.discount          = climada_xlsread(0,entity_temp_xls,'discount',1);

            fprintf('saving entity as %s\n',entity.assets.filename)
            entity_files{counter} = entity.assets.filename;
            save(entity.assets.filename,'entity')
        end
        
        if add_measures
            fprintf('adding measures\n')
            load(entity_file_mat)
            load([entities_dir filesep '20151116_measure_project_package' filesep 'measures_project_package.mat'])
            entity.measures = measures;
            entity.assets.filename = entity_file_mat;
            save(entity.assets.filename,'entity')
        end
        
    end %sheets
end % scenario_eco_i
clear entity_file_xls sheets ul_loc s_i tc_ndx damfun_file_xls entity_temp_xls fN yr_
clear fld_i flds nan_ndx entity_file_mat force_damfun_re_read force_assets_re_read asci_file



%% create measures to represent project package (20151117)
% create four measures
%  1. reduced FL depth due to a package of measures
%  2. regional scope of land raising on top of reduced FL depth, read regional scope from excel file, given from Ecorys  
%  3. reduced FL duration due to a package of measures
%  4. TC with no measures so that it runs through the calculation nontheless
clear measures
n_measures = 4;
measures = climada_measures_construct('',n_measures);

% reduced FL depth hazard due to a package of measures
asci_path = [barisal_data_dir filesep 'entities' filesep '20151116_measure_project_package' filesep 'reduced_FL_hazard'];
measure_i = 1; 
measures.name{measure_i} = 'FL depth project package';
measures.peril_ID(measure_i) = {'FL'};
measures.hazard_event_set{measure_i} = [asci_path filesep 'CaseDifferenceDepth_m.asc'];
measures.hazard_event_set_operator{measure_i} = 'plus';

% land raising at selected locations
measure_i = 2; 
measures.name{measure_i} = 'land raising';
measures.hazard_event_set{measure_i} = [asci_path filesep 'CaseDifferenceDepth_m.asc'];
measures.hazard_event_set_operator{measure_i} = 'plus';
measures.hazard_intensity_impact_b(measure_i) = -0.3; % land raising by 0.3m
measures.peril_ID(measure_i) = {'FL'};
measure_location = climada_xlsread('no',...
    [entities_dir filesep '20151116_measure_project_package' filesep 'Locations for resilient buildings measure.xlsx'],'Assets',1);
n_categories = length(entity.assets.lon)/length(measure_location.flood_resilient_buildings);
% replicate matrix for all categories
regional_scope = repmat(measure_location.flood_resilient_buildings,n_categories,1);
fprintf('Measure flood resilient buildlings\n')
fprintf('\t - %d locations (1 category), %d locations with land raising \n',length(measure_location.flood_resilient_buildings),sum(measure_location.flood_resilient_buildings))
fprintf('\t - %d all locations (%d categories) \n',length(entity.assets.lon),n_categories)
% initialize logical index to define the regional scope of measures
measures.regional_scope = ones(length(entity.assets.Value),n_measures);
% create logical 
measures.regional_scope(:,measure_i) = logical(regional_scope);

% reduced FL duration hazard due to a package of measures
measure_i = 3; 
measures.name{measure_i} = 'FL duration project package';
measures.peril_ID(measure_i) = {'FL_'};
measures.hazard_event_set{measure_i} = [asci_path filesep 'CaseDifferenceDuration_fraction.asc'];
measures.hazard_event_set_operator{measure_i} = 'times';

% TC with no measures so that it runs through the calculation nontheless
measure_i = 4; 
measures.name{measure_i} = 'TC control';
measures.peril_ID(measure_i) = {'TC'};

save([entities_dir filesep '20151116_measure_project_package' filesep 'measures_project_package.mat'],'measures')

% measures = climada_measures_check(measures,entity.assets);
% entity.measures = measures;



%% measures impact calc for today, 2030 and 2050 including climate change
% calculate measures_impact for each peril in different scenarios specified below. 
%================

peril_IDs   = {'FL_depth_monsoon' 'FL_duration_monsoon' 'FL_depth_cyclone' 'FL_duration_cyclone' 'TC'};
% peril_IDs   = {'TC'};

eco_scen    = 1; %[1 2 3];
cc_scen     = 'moderate'; %{'moderate' 'extreme'}; %
year_i      = 2014;
year_f      = [2030 2050];

silent_mode = 0;
sanity_check = 1;

%init
clear measures_impact1 measures_impact3 measures_impact5

for scenario_eco_i = 1:numel(eco_scen)
    for p_i = 1:numel(peril_IDs); %peril_IDs(1) %peril_ID = peril_IDs %
        peril_ID = peril_IDs{p_i};
        if strcmp(peril_ID,'TC')
            peril = 'cyclone';
        else
            peril = 'flood';
        end
        ed_i = p_i; % counter for EDS entries

        year_f = 2030;
        climada_global.future_reference_year = year_f;

        % EDS1 for scenario hazard and entity in present reference year
        [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files); 
        [entity,e_i] = barisal_get_entity(year_i,peril,entity_files,eco_scen(scenario_eco_i));
        hazard.scenario = 'no change';
        entity.assets.region = sprintf('BCC %s', peril_ID);
        fprintf('\n***** Scenario %d, EDS1 for %s | %s *****\n', scenario_eco_i,char(entity.assets.comment),char(strtok(hazard.comment,',')))
        scen_name1 = ['Today''s; expected damage; ' num2str(year_i)];
        measures_impact1(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
        measures_impact1(ed_i).filename = [results_dir filesep 'BCC_measure_package_impact_' num2str(year_i) '.mat'];
        for i = 1:length(measures_impact1(ed_i).EDS)
            % convert back to UTM
            [measures_impact1(ed_i).EDS(i).assets.X,measures_impact1(ed_i).EDS(i).assets.Y] = ...
                ll2utm_shift(measures_impact1(ed_i).EDS(i).assets.lat,measures_impact1(ed_i).EDS(i).assets.lon);
        end
        save(measures_impact1(ed_i).filename,'measures_impact1')   
        

        % EDS3 for climate change scenario: future entity 2030, future hazard 2030
        if scenario_eco_i == 2
            cc_scen = 'extreme';
        else 
            cc_scen = 'moderate';
        end
        [hazard,h_i] = barisal_get_hazard(year_f,cc_scen,peril_ID,hazard_files);
        [entity,e_i] = barisal_get_entity(year_f,peril,entity_files,eco_scen(scenario_eco_i));
        hazard.scenario = sprintf('%s change',cc_scen);
        entity.assets.region = sprintf('BCC %s', peril_ID);
        fprintf('\n***** Scenario %d, EDS3 for %s | %s *****\n',scenario_eco_i,char(entity.assets.comment),char(strtok(hazard.comment,',')))
        scen_name3 = sprintf('Increase;from %s;climate change;%d',cc_scen,num2str(year_f));
        measures_impact3(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
        measures_impact3(ed_i).filename = [results_dir filesep 'BCC_measure_package_impact_cc_' cc_scen '_' num2str(year_f) '.mat'];
        for i = 1:length(measures_impact3(ed_i).EDS)
            % convert back to UTM
            [measures_impact3(ed_i).EDS(i).assets.X,measures_impact3(ed_i).EDS(i).assets.Y] = ...
                ll2utm_shift(measures_impact3(ed_i).EDS(i).assets.lat,measures_impact3(ed_i).EDS(i).assets.lon);
        end
        save(measures_impact3(ed_i).filename,'measures_impact3')

        year_f = 2050;
        climada_global.future_reference_year = year_f;

        % EDS5 for climate change scenario: future entity 2050, future hazard 2050
        if scenario_eco_i == 2
            cc_scen = 'extreme';
        else 
            cc_scen = 'moderate';
        end
        [hazard,h_i] = barisal_get_hazard(year_f,cc_scen,peril_ID,hazard_files);
        [entity,e_i] = barisal_get_entity(year_f,peril,entity_files,eco_scen(scenario_eco_i));
        hazard.scenario = sprintf('%s change',cc_scen);
        entity.assets.region = sprintf('BCC %s', peril_ID);
        fprintf('\n***** Scenario %d, EDS5 for %s | %s *****\n',scenario_eco_i,char(entity.assets.comment),char(strtok(hazard.comment,',')))
        %scen_name5 = ['Increase; from ' cc_scen '; climate change; ' num2str(year_f)];
        scen_name5 = sprintf('Increase;from %s;climate change;%d',cc_scen,num2str(year_f));

        measures_impact5(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
        measures_impact5(ed_i).filename = [results_dir filesep 'BCC_measure_package_impact_cc_' cc_scen '_' num2str(year_f) '.mat'];
        for i = 1:length(measures_impact5(ed_i).EDS)
            % convert back to UTM
            [measures_impact5(ed_i).EDS(i).assets.X,measures_impact5(ed_i).EDS(i).assets.Y] = ...
                ll2utm_shift(measures_impact5(ed_i).EDS(i).assets.lat,measures_impact5(ed_i).EDS(i).assets.lon);
        end
        save(measures_impact5(ed_i).filename,'measures_impact5')
    end
end %scenario_eco_i



%% add baseline results to measures_impact
load(['\\CHRB1065.CORP.GWPNET.COM\homes\X\S3BXXW\Documents\lea\climada_git\climada_data_barisal\results\20150916_baseline' filesep 'EDS_scenario_1'])
% load([results_dir filesep 'EDS_scenario_1'])

measures_impact1_baseline = climada_measures_impact_add(measures_impact1,EDS1,entity);
measures_impact3_baseline = climada_measures_impact_add(measures_impact3,EDS3,entity);
measures_impact5_baseline = climada_measures_impact_add(measures_impact5,EDS5,entity);

ed_i = 1;
save(measures_impact1_baseline(ed_i).filename,'measures_impact1_baseline')
save(measures_impact3_baseline(ed_i).filename,'measures_impact3_baseline')
save(measures_impact5_baseline(ed_i).filename,'measures_impact5_baseline')

measures_impact = measures_impact1_baseline;
measures_impact(6:10) = measures_impact3_baseline;
measures_impact(11:15) = measures_impact5_baseline;
measures_impact(1).filename = [results_dir filesep 'measures_impact_2014_2030_2050_scenario_1']
save(measures_impact(1).filename,'measures_impact')


%% combine measures_impact for all perils per scenario (timehorizon)

peril_list = {'FL monsoon' 'FL monsoon duration'};
measures_impact = climada_measures_impact_combine_scenario(measures_impact1_baseline,measures_impact3_baseline,measures_impact5_baseline,peril_list);


% measures_impact = measures_impact1_baseline;
% measures_impact(6:10) = measures_impact3_baseline;
% measures_impact(11:15) = measures_impact5_baseline;
% 
% % get all scenario names
% for s_i=1:numel(measures_impact)
%     scenario_all{s_i,1} = measures_impact(s_i).scenario.name;
% end
% scenario_unique = unique(scenario_all);
% 
% % loop over the unique scenarios
% for s_i = 1:numel(scenario_unique)
%     is_scenario = strcmp(scenario_unique{s_i},scenario_all);
%     measures_impact_temp = [];
%     measures_impact_temp = measures_impact(is_scenario);
% 
%     combine_modus = 'delete_measures';
%     measures_impact_combined(s_i) = climada_measures_impact_combine(measures_impact_temp(1),measures_impact_temp(2:end),combine_modus);
% end



%% create some benefit plots
measure_no = 1;
fieldname_to_plot = {'ED_at_centroid' 'benefit'}; plot_method= 'plotclr'; 
timehorizon = 1;% time horizon
category_criterium = '';
[~, fig] = climada_map_plot(measures_impact,fieldname_to_plot,plot_method,measure_no,timehorizon,category_criterium);


measure_no = 3;
fieldname_to_plot = {'ED_at_centroid' 'benefit'}; plot_method= 'plotclr'; 

measure_no = 4;
fieldname_to_plot = {'ED_at_centroid'};
struct_no = 1;% peril type
category_criterium = '';
[~, fig] = climada_map_plot(measures_impact1_baseline,fieldname_to_plot,plot_method,measure_no,struct_no,category_criterium);

category_criterium = 'industry_buildings_Semi_Pucca_ASSETS_30_cm_elevation_';
category_criterium = 'Residential_buildings_Pucca_ASSETS';
category_criterium = categories(3:4);
entity.assets.Category = unique(entity.assets.Category);
measure_no = 1;
fieldname_to_plot = 'Value'; plot_method= 'plotclr'; 
struct_no = 1;% peril type
[~, fig] = climada_map_plot(entity,fieldname_to_plot,plot_method,measure_no,struct_no,category_criterium);

event_no = 1721;
fieldname_to_plot = 'intensity'; plot_method= 'plotclr'; 
struct_no = 1;% peril type
[~, fig] = climada_map_plot(hazard,fieldname_to_plot,plot_method,event_no,struct_no);

event_no = 1;
fieldname_to_plot = 'intensity'; plot_method= 'plotclr'; 
struct_no = 1;% peril type
[~, fig] = climada_map_plot('Barisal_BCC_hazard_FL_depth_cyclone_2014',fieldname_to_plot,plot_method,event_no,struct_no);

event_no = 1721;
fieldname_to_plot = 'intensity'; plot_method= 'contourf'; 
struct_no = 1;% peril type
[~, fig] = climada_map_plot('Barisal_BCC_hazard_TC_2014',fieldname_to_plot,plot_method,event_no,struct_no);



%% write ED output for measures_impact

clear EDS_project_package_1 EDS_project_package_3 EDS_project_package_5 EDS_project_package

% collect EDS for all perils (FL depth, FL duration, FL depth cyclone, FL duration cyclone, TC)
EDS_project_package_1 = [measures_impact1(1).EDS measures_impact1(2).EDS measures_impact1(3).EDS...
                        measures_impact1(4).EDS measures_impact1(5).EDS(1)];
% EDS_project_package_1 = [measures_impact1(3).EDS measures_impact1(4).EDS];    
EDS_project_package_3 = [measures_impact3(1).EDS measures_impact3(2).EDS measures_impact3(3).EDS...
                        measures_impact3(4).EDS measures_impact3(5).EDS(1)];
EDS_project_package_5 = [measures_impact5(1).EDS measures_impact5(2).EDS measures_impact5(3).EDS...
                        measures_impact5(4).EDS measures_impact5(5).EDS(1)];
% convert back to UTM
for ed_i = 1:length(EDS_project_package_1)
    [EDS_project_package_1(ed_i).assets.X,EDS_project_package_1(ed_i).assets.Y] = ...
                   ll2utm_shift(EDS_project_package_1(ed_i).assets.lat, EDS_project_package_1(ed_i).assets.lon);
    [EDS_project_package_3(ed_i).assets.X,EDS_project_package_3(ed_i).assets.Y] = ...
                   ll2utm_shift(EDS_project_package_3(ed_i).assets.lat, EDS_project_package_3(ed_i).assets.lon);
    [EDS_project_package_5(ed_i).assets.X,EDS_project_package_5(ed_i).assets.Y] = ...
                   ll2utm_shift(EDS_project_package_5(ed_i).assets.lat, EDS_project_package_5(ed_i).assets.lon);           
    EDS_project_package_1(ed_i).peril_ID = '';
    EDS_project_package_3(ed_i).peril_ID = '';
    EDS_project_package_5(ed_i).peril_ID = '';
end                    
                
% loop over three time horizons
timehorizons = [2014 2030 2050];                    
for t_i = 1%:3

    % fill variable with a specific time horizon
    switch t_i
        case 1
            EDS_project_package = EDS_project_package_1;
        case 2 
            EDS_project_package = EDS_project_package_3;
        case 3
            EDS_project_package = EDS_project_package_5;
    end
    
    % write ED per category report
    EDS_report_xls = [results_dir filesep sprintf('BCC_ED_report_scenario_%d_project_package_%s.xls',eco_scen(scenario_eco_i),datestr(now,'yyyymmdd'))];
    %if exist(EDS_report_xls,'file'), delete(EDS_report_xls); end
    benefit_flag = 0;
    assets_flag = 1;
    sheetname = sprintf('ED_per_category_%d',timehorizons(t_i));
    output_report = climada_EDS_ED_per_category_report(entity, EDS_project_package, EDS_report_xls,sheetname,benefit_flag,0,assets_flag);
    
    % write ED at centroid to excel
    EDS_report_xls = [results_dir filesep sprintf('BCC_ED_report_scenario_%d_project_package_at_centroid_%s.xls',eco_scen(scenario_eco_i),datestr(now,'yyyymmdd'))];
    sheetname = sprintf('ED_at_centroid_%d',timehorizons(t_i));
    climada_EDS_ED_at_centroid_report_xls(EDS_project_package,EDS_report_xls,sheetname);
end





%%
figure
assets_indx = strcmp(entity.assets.Category,'Residential_buildings_Pucca_ASSETS');
plotclr(entity.assets.lon(assets_indx), entity.assets.lat(assets_indx), AED_rel(assets_indx),'','',1,0.01,5)
title('AED relative, Residential_buildings_Pucca_ASSETS')

figure
assets_indx = strcmp(entity.assets.Category,'Residential_buildings_Pucca_ASSETS');
plotclr(entity.assets.lon(assets_indx), entity.assets.lat(assets_indx), AED_rel(assets_indx),'','',1)
title('AED relative, Residential_buildings_Pucca_ASSETS')


figure
assets_indx = strcmp(entity.assets.Category,'Residential_buildings_Katcha_ASSETS');
plotclr(entity.assets.lon(assets_indx), entity.assets.lat(assets_indx), AED_rel(assets_indx),'','',1,0.01,5)
title('AED relative, Residential_buildings_Katcha_ASSETS')


figure
assets_indx = strcmp(entity.assets.Category,'Residential_buildings_Pucca_ASSETS');
plotclr(entity.assets.lon(assets_indx), entity.assets.lat(assets_indx), entity.assets.income(assets_indx),'','',1)
title('Income, Residential_buildings_Pucca_ASSETS')

figure
assets_indx = strcmp(entity.assets.Category,'Residential_buildings_Pucca_ASSETS');
plotclr(entity.assets.lon(assets_indx), entity.assets.lat(assets_indx), measures_impact3(1).EDS(1).ED_at_centroid(assets_indx),'','',1)
title('ED, Residential_buildings_Pucca_ASSETS')

return


%% measures construction

% structure measures, assign measures.hazard_event_set by searching through
% barisal entities folder and finding measure folders starting with the
% substring "Measures"

e_dir_ = dir(entities_dir);

%init
clear measures; measures = climada_measures_construct([],0);

% hazard modifying measures
for e_dir_i = 1:length(e_dir_)
    if e_dir_(e_dir_i).isdir && ...
            ~isempty(strfind(upper(e_dir_(e_dir_i).name),'MEASURES'))
        
        e_subdir_ = dir([entities_dir filesep e_dir_(e_dir_i).name]);
        
        for m_i = 1:length(e_subdir_)
            if e_subdir_(m_i).name(1) == '.', continue; end
            asci_file = [entities_dir filesep e_dir_(e_dir_i).name filesep e_subdir_(m_i).name];
            % construct name from folder and file
            name    = strrep([e_dir_(e_dir_i).name],'Measures_','');
            name(1) = upper(name(1)); name = strrep(name,'_',' ');
            
            measures = climada_measures_construct(measures,1);
            
            measures.name{end}              = name;
            measures.hazard_event_set{end}  =  asci_file;
            % operator required for call to climada_distributed_measures in
            % climada_measures_impact_advanced. Make sure it is correct by
            % checking peril ID => multiplication for flood duration and
            % addition for flood depth
            if ~isempty(strfind(e_subdir_(m_i).name,'Duration'))
                measures.peril_ID{end} = 'FL_';
                measures.hazard_event_set_operator{end} = 'times';
            else
                measures.peril_ID{end} = 'FL';
                measures.hazard_event_set_operator{end} = 'plus';
            end
            % fprintf('%s \t%s\n',name,measures.peril_ID{end})
        end
    end
end

% % assign measures struct to entities and save each one again
measures                              = climada_measures_construct(measures,4);
measures.name{end-3}                  = 'Flood resilient crops';
measures.damagefunctions_map{end-3}   = '1320to1321';
measures.peril_ID{end-3}              = 'FL_';

% save measures struct (before adding entity file)
measures.filename = [entities_dir filesep 'BCC_measures_' datestr(now,'ddmmyy') '.mat'];
fprintf('saving measures as %s\n',measures.filename)
save(measures.filename,'measures')

measures_ori = measures; % for ponds retain routine

for e_i = 1:length(entity_files)
    
    measures = measures_ori;
    
    load(entity_files{e_i})
    
    % remove old measures, just in case
    if isfield(entity,'measures'),entity = rmfield(entity,'measures'); end
    
    measures.name{end-2}        = 'Spatial planning';
    measures.entity_file{end-2} = [fileparts(entity_files{e_i}) filesep...
        'Measure_Spatial_planning_290615_se_1_' lower(strrep(entity.assets.comment,' ','_')) '.mat'];

    if ~exist(measures.entity_file{end-2},'file')
        cprintf([ 1 0.5 0],'WARNING: entity file for measure %s does not exist\n',measures.name{end-2})
    end
    
    measures.name{end-1}        = 'Early warning system';
    measures.entity_file{end-1} = [fileparts(entity_files{e_i}) filesep...
        'Measure_Early_warning_system_290615_se_1_' lower(strrep(entity.assets.comment,' ','_')) '.mat'];
    measures.damagefunctions_map{end-1}   = '1210to1211;1220to1221';

    if ~exist(measures.entity_file{end-1},'file')
        cprintf([ 1 0.5 0],'WARNING: entity file for measure %s does not exist\n',measures.name{end-1})
    end
    measures.name{end}          = 'Flood proof road infrastructure';
    measures.entity_file{end}   = [fileparts(entity_files{e_i}) filesep...
        'Measure_Flood_proof_road_infrastructure_290615_se_1_' lower(strrep(entity.assets.comment,' ','_')) '.mat'];
 
    if ~exist(measures.entity_file{end},'file')
        cprintf([ 1 0.5 0],'WARNING: entity file for measure %s does not exist\n',measures.name{end})
    end   
    % drainage for correct refence year
    rm_ndx = [];
    for m_i = 1:length(measures.name)
        if ~isempty(strfind(measures.name{m_i},'Ponds retain'))
            if isempty(strfind(measures.name{m_i},num2str(entity.assets.reference_year)))
                rm_ndx = [rm_ndx m_i];
            end
            measures.name{m_i} = 'Ponds retain';
        end
    end
    
    entity.measures = climada_measures_construct(measures,-rm_ndx);
%     entity.measures = measures;
    fprintf('measures assigned to entity %s:\n',entity.assets.comment)
    for measure_i = 1:length(entity.measures.name)
        cprintf(entity.measures.color_RGB(measure_i,:),'\t%s\n',entity.measures.name{measure_i});
    end
    
    fprintf('saving entity with measures to %s\n',entity.assets.filename)
    
    save(entity.assets.filename,'entity')
end
% display measures in their respective colours
% for measure_i = 1:length(measures.name),cprintf(measures.color_RGB(measure_i,:),'%s\n',measures.name{measure_i}); end
clear e_i e_dir_i e_dir_ e_subdir_ e_subdir_ m_i m_i measure_i rm_ndx name

%% measures package construction

%init
clear measures; measures = climada_measures_construct([],3);

% hazard modifying measures

asci_path = [barisal_data_dir filesep 'entities' filesep 'Measures_package'];

% flood depth
measures.name{1}                =   'Measure package';
measures.hazard_event_set{1}    =   [asci_path filesep 'CaseDifferenceDepth_m.asc'];
measures.peril_ID{1}            =   'FL';
measures.hazard_event_set_operator{1} = 'plus';
measures.damagefunctions_map{1}   = '1210to1211;1220to1221'; % vehicles

% flood duration
measures.name{2}                =   'Measure package';
measures.hazard_event_set{2}    =   [asci_path filesep 'CaseDifferenceDuration_fraction.asc'];
measures.peril_ID{2}            =   'FL_';
measures.hazard_event_set_operator{2} = 'times';
measures.damagefunctions_map{2} = '1320to1321'; % Flood resilient crops (only for duration)

% cyclone wind
measures.name{3}                =   'Measure package';
measures.peril_ID{3}            =   'TC';

% save measures struct (before adding entity file)
measures.filename = [entities_dir filesep 'BCC_measure_package' datestr(now,'ddmmyy') '.mat'];
fprintf('saving measures as %s\n',measures.filename)
save(measures.filename,'measures')

for y_i = [2014 2030 2050]

    % set correct measure entity file
    measures.entity_file{1} = [entities_dir filesep...
        'Measure_package_070715_se_1_floods_' num2str(y_i) '.mat'];
    measures.entity_file{2} = [entities_dir filesep...
        'Measure_package_070715_se_1_floods_' num2str(y_i) '.mat'];
    measures.entity_file{3} = [entities_dir filesep...
        'Measure_package_070715_se_1_cyclones_' num2str(y_i) '.mat'];
    
    % assign measures to baseline entity struct
    % for both flood depth and duration
    [entity, e_i] = barisal_get_entity(y_i,'floods',entity_files);
    % remove old measures, just in case
    if isfield(entity,'measures'),entity = rmfield(entity,'measures'); end
    entity.measures = measures;
    fprintf('saving entity with measures to %s\n',entity.assets.filename)
    save(entity.assets.filename,'entity')
    
    % for TC
    [entity, e_i] = barisal_get_entity(y_i,'cyclones',entity_files);
    % remove old measures, just in case
    if isfield(entity,'measures'),entity = rmfield(entity,'measures'); end
    entity.measures = measures;
    fprintf('saving entity with measures to %s\n',entity.assets.filename)
    save(entity.assets.filename,'entity')
end
% display measures in their respective colours
for measure_i = 1:length(measures.name),cprintf(measures.color_RGB(measure_i,:),'%s\n',measures.name{measure_i}); end
clear e_i e_dir_i e_dir_ e_subdir_ e_subdir_ m_i m_i measure_i rm_ndx name



%% population entity & casualties [uncomment, but be careful]
% 
% Entity excel file from Ecorys (unstructured)
entity_file_xls = [entities_dir filesep 'BCC_entity_210615_scenario 1_population.xls'];

% Entity template file from global data dir
entity_temp_xls = [climada_global.data_dir filesep 'entities' filesep 'entity_template.xls'];
sheets          = {'2014' '2030' '2050'}; % for population entity

% load damagefunctions
load([entities_dir filesep 'BCC_population_entity_090615.mat'])
damagefunctions = entity.damagefunctions;

force_assets_re_read   = 1;% Whether to re-read entity from xls (1) or load matfiles (0)

for s_i = 1:length(sheets)
    clear entity
    
    % mat file name - separate mat file for each sheet in entity xls
    entity_file_mat     = strrep(entity_file_xls,'.xls',['_' lower(sheets{s_i}) '.mat']);
    
    if exist(entity_file_mat,'file') && ~force_assets_re_read
        % Load mat file
        fprintf('entity %s .mat file already exists, skipping\n',lower(sheets{s_i}))
        entity_files{s_i} = entity_file_mat;
        load(entity_file_mat)
        % append filename for consistency
        entity.assets.filename = entity_file_mat;
        % save corrected entity
        save(entity.assets.filename,'entity')
    else
        % read in entity from xls
        entity = climada_entity_read(entity_temp_xls);
        %entity = climada_entity_read_wo_assets(entity_temp_xls);
        entity.assets = climada_xlsread('no',entity_file_xls,sheets{s_i},1);
        entity.assets.comment = ['Population entity ' sheets{s_i}];
        entity.assets.reference_year = str2num(sheets{s_i});
        % coord transformation from UTM to lat lon
        [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.lon, entity.assets.lat);
        entity.assets.filename          = entity_file_mat; % for filename consistency
        % use people damagefunctions
        entity.damagefunctions = damagefunctions;
        fprintf('saving entity as %s\n',entity.assets.filename)
        save(entity.assets.filename,'entity')
        entity_files{s_i} = entity.assets.filename;
    end
    
end
clear entity_file_xls sheets s_i entity_temp_xls entity_file_mat force_assets_re_read


%%
climada_global.Value_unit = 'people';

year_i = 2014; climada_global.present_reference_year = year_i;
year_f = 2030; climada_global.future_reference_year  = year_f;
cc_scen = 'moderate'; ed_i = 0; clear EDS1 EDS2 EDS3 % init

for peril_ID = {'FL_depth_monsoon' 'FL_depth_cyclone'}
    ed_i = ed_i +1;
    % EDS1 for scenario hazard and entity in present reference year
    [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
    [entity,e_i] = barisal_get_entity(year_i,'population',entity_files);
    fprintf('***** EDS1 for %s | %s *****\n',char(entity.assets.comment),char(strtok(hazard.comment,',')))
    EDS1(ed_i)    = climada_EDS_calc(entity,hazard,'',1,'',1);
    scen_name1 = ['Today''s; expected casualties; ' num2str(year_i)];
    fprintf('Annual expected casualties: %2.2f mn\n',EDS1(ed_i).ED)
    
    % EDS2 economic growth, no climate change
    [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
    [entity,e_i] = barisal_get_entity(year_f,'population',entity_files);
    fprintf('***** EDS2 for %s | %s *****\n',char(entity.assets.comment),char(strtok(hazard.comment,',')))
    EDS2(ed_i)    = climada_EDS_calc(entity,hazard,'',1);
    scen_name2 = ['Increase; from population; growth ' num2str(year_f)];
    fprintf('Annual expected casualties: %2.2f mn\n',EDS2(ed_i).ED)
    
    % EDS3 climate change + economic growth
    [hazard,h_i] = barisal_get_hazard(year_f,cc_scen,peril_ID,hazard_files);
    [entity,e_i] = barisal_get_entity(year_f,'population',entity_files);
    fprintf('***** EDS3 for %s | %s *****\n',char(entity.assets.comment),char(strtok(hazard.comment,',')))
    EDS3(ed_i)    = climada_EDS_calc(entity,hazard,'',1);
    scen_name3 = ['Increase; from ' cc_scen '; climate change; ' num2str(year_f)];
    fprintf('Annual expected casualties: %2.2f mn\n',EDS3(ed_i).ED)
end
% figure; climada_ED_plot(EDS1, 0,'people',0,0)
% figure; climada_ED_plot(EDS2, 0,'people',0,0)
% figure; climada_ED_plot(EDS3, 0,'people',0,0)

climada_waterfall_graph_multi_peril(0,'people',EDS1,scen_name1,EDS2,scen_name2,EDS3,scen_name3)
% print(gcf,'-dpng',[results_dir filesep 'BCC_Barisal_population_waterfall_2030.png'])
% 
% 
save([results_dir filesep 'Barisal_EDS_people_2014_2030_FL_monsoon_cyclone.mat'],'EDS1', 'EDS2', 'EDS3')

% EDS(1) = climada_EDS_combine(EDS1(1),EDS1(2));

EDS(1) = EDS1(1);
EDS(1).ED = sum([EDS1.ED]);
EDS(2) = EDS2(1);
EDS(2).ED = sum([EDS2.ED]);
EDS(3) = EDS3(1);
EDS(3).ED = sum([EDS3.ED]);
EDS(1).hazard.comment = 'FL (monsoon and cyclone)';
EDS(2).hazard.comment = 'FL (monsoon and cyclone)';
EDS(3).hazard.comment = 'FL (monsoon and cyclone)';
EDS(1).assets.filename = 'BCC population 2014';
EDS(2).assets.filename = 'BCC population 2030';
EDS(3).assets.filename = 'BCC population 2030';

climada_waterfall_graph(EDS(1), EDS(2), EDS(3),'AED');
save([results_dir filesep 'Barisal_EDS_people_2014_2030.mat'],'EDS')
savefig([results_dir filesep 'Barisal_EDS_people_2014_2030.fig'])




%% Multi scenario EDS calc for waterfall

%================
% This section is key for measure impact calculations
%================

% calculate EDS for each peril in different scenarios specified below. Each
% EDSn corresponds to one bar in the waterfall graph. They should have the
% same length, and each entry in each EDSn corresponds to a peril_ID
% measure impact is also calculated, unless EDS_only is set to 1

peril_IDs   = {'FL_depth_monsoon' 'FL_duration_monsoon' 'FL_depth_cyclone' 'FL_duration_cyclone' 'TC'};
% peril_IDs   = {'TC'};

eco_scen    = 1; %[1 2 3];
cc_scen     = 'moderate'; %{'moderate' 'extreme'}; %
year_i      = 2014;
year_f      = [2030 2050];

EDS_only    = 0; %1; %EDS_only    = 0;
silent_mode = 0;
sanity_check = 1;

%init
% clear EDS1 EDS2 EDS3 EDS4 EDS5 
% clear measures_impact1 measures_impact3 measures_impact5

for scenario_eco_i = 1:numel(eco_scen)
    % for year_f = [2050 2030]
    %ed_i = 0;
    for p_i = 4%1:numel(peril_IDs); %peril_IDs(1) %peril_ID = peril_IDs %
        peril_ID = peril_IDs{p_i};
        if strcmp(peril_ID,'TC')
            peril = 'cyclone';
        else
            peril = 'flood';
        end
        %ed_i = ed_i+1; % counter for EDS entries
        ed_i = p_i; % counter for EDS entries

        year_f = 2030;
        climada_global.future_reference_year = year_f;

        % EDS1 for scenario hazard and entity in present reference year
        [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
        [entity,e_i] = barisal_get_entity(year_i,peril,entity_files,eco_scen(scenario_eco_i));
        %EDS1    = barisal_get_EDS(EDS,entity_files{e_i},hazard_files{h_i});
        fprintf('\n***** Scenario %d, EDS1 for %s | %s *****\n', scenario_eco_i,char(entity.assets.comment),char(strtok(hazard.comment,',')))
        scen_name1 = ['Today''s; expected damage; ' num2str(year_i)];
        %fprintf('Annual expected damage: %2.2f mn\n',EDS1(ed_i).ED/1000000)
        if ~EDS_only
            measures_impact1(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
            measures_impact1(ed_i).filename = [results_dir filesep 'BCC_measure_package_impact_' num2str(year_i) '.mat'];
            for i = 1:length(measures_impact1(ed_i).EDS)
                % convert back to UTM
                [measures_impact1(ed_i).EDS(i).assets.X,measures_impact1(ed_i).EDS(i).assets.Y] = ...
                    ll2utm_shift(measures_impact1(ed_i).EDS(i).assets.lat,measures_impact1(ed_i).EDS(i).assets.lon);
            end
            save(measures_impact1(ed_i).filename,'measures_impact1')
        else
            EDS1(ed_i)    = climada_EDS_calc(entity,hazard,'',1,silent_mode,sanity_check);
        end        

        if EDS_only
            % EDS2 for socio-economic growth scenario: present hazard, future entity 2030
            [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
            [entity,e_i] = barisal_get_entity(year_f,peril,entity_files,eco_scen(scenario_eco_i));
            %EDS2    = barisal_get_EDS(EDS,entity_files{e_i},hazard_files{h_i});
            fprintf('\n***** Scenario %d, EDS2 for %s | %s *****\n',scenario_eco_i,char(entity.assets.comment),char(strtok(hazard.comment,',')))
            EDS2(ed_i)    = climada_EDS_calc(entity,hazard,'',1,silent_mode,sanity_check);
            scen_name2 = sprintf('Increase;from economic;growth %d;scenario %d',num2str(year_f),eco_scen(scenario_eco_i));
            fprintf('Annual expected damage: %2.2f mn\n',EDS2(ed_i).ED/1000000)
            %entity.measures = climada_measures_construct([],1);
            %entity.measures.name{1} = 'control';
            %measures_impact2(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
            %for i = 1:length(measures_impact2(ed_i).EDS)
            %    % convert back to UTM
            %    [measures_impact2(ed_i).EDS(i).assets.X,measures_impact2(ed_i).EDS(i).assets.Y] = ...
            %        ll2utm_shift(measures_impact2(ed_i).EDS(i).assets.lat,measures_impact2(ed_i).EDS(i).assets.lon);
            %end
        end

        % EDS3 for climate change scenario: future entity 2030, future hazard 2030
        if scenario_eco_i == 2
            cc_scen = 'extreme';
        else 
            cc_scen = 'moderate';
        end
        [hazard,h_i] = barisal_get_hazard(year_f,cc_scen,peril_ID,hazard_files);
        [entity,e_i] = barisal_get_entity(year_f,peril,entity_files,eco_scen(scenario_eco_i));
        fprintf('\n***** Scenario %d, EDS3 for %s | %s *****\n',scenario_eco_i,char(entity.assets.comment),char(strtok(hazard.comment,',')))
        %scen_name3 = ['Increase; from ' cc_scen '; climate change; ' num2str(year_f)];
        scen_name3 = sprintf('Increase;from %s;climate change;%d',cc_scen,num2str(year_f));
        %fprintf('Annual expected damage: %2.2f mn\n',EDS3(ed_i).ED/1000000)
        if ~EDS_only
            measures_impact3(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
            measures_impact3(ed_i).filename = [results_dir filesep 'BCC_measure_package_impact_cc_' cc_scen '_' num2str(year_f) '.mat'];
            for i = 1:length(measures_impact3(ed_i).EDS)
                % convert back to UTM
                [measures_impact3(ed_i).EDS(i).assets.X,measures_impact3(ed_i).EDS(i).assets.Y] = ...
                    ll2utm_shift(measures_impact3(ed_i).EDS(i).assets.lat,measures_impact3(ed_i).EDS(i).assets.lon);
            end
            save(measures_impact3(ed_i).filename,'measures_impact3')
        else
            EDS3(ed_i)    = climada_EDS_calc(entity,hazard,'',1,silent_mode,sanity_check);
        end

        year_f = 2050;
        climada_global.future_reference_year = year_f;

        if EDS_only
            % EDS4 for socio-economic growth scenario: present hazard, future entity 2050
            [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
            [entity,e_i] = barisal_get_entity(year_f,peril,entity_files,eco_scen(scenario_eco_i));
            fprintf('\n***** Scenario %d, EDS4 for %s | %s *****\n',scenario_eco_i,char(entity.assets.comment),char(strtok(hazard.comment,',')))
            EDS4(ed_i)    = climada_EDS_calc(entity,hazard,'',1,silent_mode,sanity_check);
            %scen_name4 = ['Increase; from economic; growth ' num2str(year_f)];
            scen_name4 = sprintf('Increase;from economic; growth %d;scenario %d',num2str(year_f),eco_scen(scenario_eco_i));
            fprintf('Annual expected damage: %2.2f mn\n',EDS4(ed_i).ED/1000000)
            %entity.measures = climada_measures_construct([],1);
            %entity.measures.name{1} = 'control';
            %measures_impact4(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
            %for i = 1:length(measures_impact4(ed_i).EDS)
            %    % convert back to UTM
            %    [measures_impact4(ed_i).EDS(i).assets.X,measures_impact4(ed_i).EDS(i).assets.Y] = ...
            %        ll2utm_shift(measures_impact4(ed_i).EDS(i).assets.lat,measures_impact4(ed_i).EDS(i).assets.lon);
            %end

        end

        % EDS5 for climate change scenario: future entity 2050, future hazard 2050
        if scenario_eco_i == 2
            cc_scen = 'extreme';
        else 
            cc_scen = 'moderate';
        end
        [hazard,h_i] = barisal_get_hazard(year_f,cc_scen,peril_ID,hazard_files);
        [entity,e_i] = barisal_get_entity(year_f,peril,entity_files,eco_scen(scenario_eco_i));
        fprintf('\n***** Scenario %d, EDS5 for %s | %s *****\n',scenario_eco_i,char(entity.assets.comment),char(strtok(hazard.comment,',')))
        %scen_name5 = ['Increase; from ' cc_scen '; climate change; ' num2str(year_f)];
        scen_name5 = sprintf('Increase;from %s;climate change;%d',cc_scen,num2str(year_f));
        if ~EDS_only
            measures_impact5(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
            measures_impact5(ed_i).filename = [results_dir filesep 'BCC_measure_package_impact_cc_' cc_scen '_' num2str(year_f) '.mat'];
            for i = 1:length(measures_impact5(ed_i).EDS)
                % convert back to UTM
                [measures_impact5(ed_i).EDS(i).assets.X,measures_impact5(ed_i).EDS(i).assets.Y] = ...
                    ll2utm_shift(measures_impact5(ed_i).EDS(i).assets.lat,measures_impact5(ed_i).EDS(i).assets.lon);
            end
            save(measures_impact5(ed_i).filename,'measures_impact5')
        else
            EDS5(ed_i)    = climada_EDS_calc(entity,hazard,'',1,silent_mode,sanity_check);
            fprintf('Annual expected damage: %2.2f mn\n',EDS5(ed_i).ED/1000000)
        end
    end
    %save([results_dir filesep sprintf('EDS_scenario_%d.mat',scenario_eco_i)], 'EDS1', 'EDS2', 'EDS3', 'EDS4', 'EDS5')
    clear EDS1 EDS2 EDS3 EDS4 EDS5
end %scenario_eco_i



%%
MI_EDS_report_xls = [results_dir filesep 'BCC_ED_report_Measure_package_' datestr(now,'ddmmyy') '.xls'];
if ~EDS_only
    for m = [1 3 5]
        measures_impact = eval(['measures_impact' num2str(m)]);
%         eval(['MI_EDS_combined' num2str(m) ' = climada_measures_impact_report(measures_impact,''NO_SAVE''); '])
        eval(['MI_EDS_combined' num2str(m) ' = climada_measures_impact_report(measures_impact,MI_EDS_report_xls); '])
    end
end

%%

% plotting frenzy

% for ed_i = 1:length(EDS1)
%     climada_ED_plot(EDS1(ed_i), 0,'BDT',0,0)
%     print(gcf,'-dpng',[results_dir filesep 'BCC_dmg_slums_' strrep(strtok(EDS1(ed_i).hazard.comment,','),' ','_') '.png'])
%     shape_plotter(BCC_wards_ll,'','linewidth',1,'color',[81 81 81]/255);
%     close
%     climada_ED_plot(EDS3(ed_i), 0,'BDT',0,0)
%     print(gcf,'-dpng',[results_dir filesep 'BCC_dmg_slums_' strrep(strtok(EDS3(ed_i).hazard.comment,','),' ','_') '.png'])
%     shape_plotter(BCC_wards_ll,'','linewidth',1,'color',[81 81 81]/255);
%     close
%     climada_ED_plot(EDS5(ed_i), 0,'BDT',0,0)
%     print(gcf,'-dpng',[results_dir filesep 'BCC_dmg_slums_' strrep(strtok(EDS5(ed_i).hazard.comment,','),' ','_') '.png'])
%     shape_plotter(BCC_wards_ll,'','linewidth',1,'color',[81 81 81]/255);
%     close
% end


%% measures impact (benefits) report per peril
barisal_MI_per_peril(measures_impact5,measures,peril_IDs)


%% multi peril waterfall graph
year_i = 2014;
scen_name1 = ['Today''s; expected damage; ' num2str(year_i)];
year_f = 2030;
eco_scen    = [1 2 3];scenario_eco_i = 2;
scen_name2 = sprintf('Increase;from economic;growth %d;scenario %d',num2str(year_f),eco_scen(scenario_eco_i));
cc_scen = 'moderate';
scen_name3 = sprintf('Increase;from %s;climate change;%d',cc_scen,num2str(year_f));

if exist('EDS1','var') && exist('EDS2','var') && exist('EDS3','var') && exist('EDS4','var') && exist('EDS5','var')
    % multi peril waterfall 2030
    fig = climada_waterfall_graph_multi_peril(0,'BDT',EDS1,scen_name1,EDS2,scen_name2,EDS3,scen_name3);
    print(fig,'-dpng',[results_dir filesep 'BCC_waterfall_multi_peril_2030_' char(cc_scen) '.png'])

    % multi peril waterfall 2050
    fig = climada_waterfall_graph_multi_peril(0,'BDT',EDS1,scen_name1,EDS4,scen_name4,EDS5,scen_name5);
    print(fig,'-dpng',[results_dir filesep 'BCC_waterfall_multi_peril_2050_' char(cc_scen) '.png'])
    
    % multi horizon waterfall (watch out! may be a bit dodgy...)
    %fig = climada_waterfall_graph_2timehorizons('AED',0,'BDT',...
    %    MI_EDS_combined1(end),   ...
    %    MI_EDS_combined2(end),   ...
    %    MI_EDS_combined3(end),   ...
    %    MI_EDS_combined4(end),   ...
    %    MI_EDS_combined5(end)          );
    %print(fig,'-dpng',[results_dir filesep 'BCC_waterfall_multi_horizon_' char(cc_scen) '.png'])
end

for s_i = [1 3 5]
    MI_EDS_combined = eval(['MI_EDS_combined' num2str(s_i)]);
    for ed_i = 1:length(MI_EDS_combined)
        % use UTM X/Y instead of lat/lon, temporarily overwrite lat/lon
        MI_EDS_combined(ed_i).assets.lon = MI_EDS_combined(ed_i).assets.X;
        MI_EDS_combined(ed_i).assets.lat = MI_EDS_combined(ed_i).assets.Y;
        climada_MI_plot(MI_EDS_combined(ed_i), 0,'BDT',0,0,1)
        shape_plotter(BCC_wards_ll,'','x','y','linewidth',1,'color',[81 81 81]/255);
        print(gcf,'-dpng',[results_dir filesep 'BCC_measure_benefit_' ...
            strrep(MI_EDS_combined(ed_i).annotation_name,' ','_') '_' ...
            num2str(MI_EDS_combined(1).reference_year) '.png'])
        close
    end
end

clear -regexp ^scen_name\d{1}$
clear m measures_impact cc_scen ed_i s_i 

%% baseline EDS for ED report for Ecorys
% EDS_baseline = [EDS1 EDS2 EDS3 EDS4 EDS5]; save([results_dir filesep 'EDS_baseline.mat'],'EDS_baseline')

climada_global.present_reference_year = 2014;

eco_scen    = [1 2 3];
cc_scen     = 'moderate';

for scenario_eco_i = 1:numel(eco_scen)
    clear EDS_baseline 
    
    load([results_dir filesep sprintf('EDS_scenario_%d.mat',eco_scen(scenario_eco_i))])
    
    if scenario_eco_i == 2
        cc_scen = 'extreme';
    else 
        cc_scen = 'moderate';
    end
    % annotation 
    year_i = 2014;
    scen_name1 = ['Today''s; expected damage; ' num2str(year_i)];
    year_f = 2030;
    scen_name2 = sprintf('Increase\nfrom economic\ngrowth %d,\nscenario %d',year_f,eco_scen(scenario_eco_i));
    scen_name3 = sprintf('Increase\nfrom %s\nclimate change\n%d',cc_scen,year_f);
    year_f = 2050;
    scen_name4 = sprintf('Increase\nfrom economic\ngrowth %d,\nscenario %d',year_f,eco_scen(scenario_eco_i));
    scen_name5 = sprintf('Increase\nfrom %s\nclimate change\n%d',cc_scen,year_f);
    
    % multi peril waterfall 2030
    fig = climada_waterfall_graph_multi_peril(0,'BDT',EDS1,scen_name1,EDS2,scen_name2,EDS3,scen_name3);
    print(fig,'-dpdf',[results_dir filesep sprintf('BCC_waterfall_multi_peril_scenario_%d_2030_%s.pdf',eco_scen(scenario_eco_i),char(cc_scen))])

    % multi peril waterfall 2050
    fig = climada_waterfall_graph_multi_peril(0,'BDT',EDS1,scen_name1,EDS4,scen_name4,EDS5,scen_name5);
    print(fig,'-dpdf',[results_dir filesep sprintf('BCC_waterfall_multi_peril_scenario_%d_2050_%s.pdf',eco_scen(scenario_eco_i),char(cc_scen))])

%     EDS_baseline = [EDS1 EDS3 EDS5]; %save([results_dir filesep 'EDS_baseline.mat'],'EDS_baseline')
%     EDS_report_xls = [results_dir filesep sprintf('BCC_ED_report_scenario_%d_%s.xls',eco_scen(scenario_eco_i),datestr(now,'yyyymmdd'))];
%     EDS_report_csv = [results_dir filesep sprintf('BCC_ED_report_scenario_%d_%s.csv',eco_scen(scenario_eco_i),datestr(now,'yyyymmdd'))];
% 
%     if exist(EDS_report_xls,'file'), delete(EDS_report_xls); end
%     % convert back to UTM
%     for ed_i = 1:length(EDS_baseline)
%         [EDS_baseline(ed_i).assets.X,EDS_baseline(ed_i).assets.Y] = ...
%             ll2utm_shift(EDS_baseline(ed_i).assets.lat,EDS_baseline(ed_i).assets.lon);
%         EDS_baseline(ed_i).peril_ID = '';
%     end
% 
%     % write to csv
%     % climada_EDS_ED_at_centroid_report(EDS_baseline,EDS_report_csv);
% 
%     % write to excel
%     climada_EDS_ED_at_centroid_report_xls(EDS_baseline,EDS_report_xls,'ED_at_centroid');
%     
%     % write ED per category report
%     benefit_flag = 0;
%     assets_flag = 1;
%     output_report = climada_EDS_ED_per_category_report(entity, EDS_baseline, EDS_report_xls,'ED_per_category',benefit_flag,0,assets_flag);
end

clear peril_ID peril year_i year_f EDS_only EDS_report_xls h_i e_i peril_IDs CC_SCEN YEAR_F EDS_report_csv

%% damage calc

% This section has become somewhat redundant.

% calculate EDS for each entity-hazard pair, where they have coinciding
% reference years and the correct peril_ID matching. Store as struct array
% and save as EDS_save_file
EDS_save_file = [results_dir filesep 'BCC_EDS_' datestr(now,'ddmmyy') '.mat'];
EDS_load_file = [results_dir filesep 'BCC_EDS_090615.mat'];

EDS_force_recalc = 0;

if exist(EDS_load_file,'file') && ~EDS_force_recalc
    load(EDS_load_file)
else
    EDS        = climada_EDS_multi_calc(entity_files,hazard_files,EDS_save_file,1,0);
end

% But this may still be useful to get the total asset base into the EDS
% struct for climada_ED_plot
% total value from max of flood entity and cyclone entity at each point
for i = 1:5
    EDS = eval(['EDS' num2str(i)]);
    for ed_i = 1:length(EDS)
        ndx = find([EDS.reference_year] == EDS(ed_i).reference_year);
        max_val = [];
        for ndx_i = ndx
            max_val = max([max_val EDS(ndx_i).assets.Value],[],2);
        end
        
        EDS(ed_i).Value_total = sum(max_val);
    end
    eval(['EDS' num2str(i) '=EDS;']);
end
clear max_val ndx EDS_force_recalc ed_i EDS_load_file EDS_save_file ndx_i


%% multi peril adaptation cost curve

% set parameters to key scenario 2030, moderate climate change, so far
% measures for flood only (TC wind excluded)
hazard_ref_year = 2030;
entity_ref_year = 2030;
cc_scen         = 'extreme'; %{'moderate' 'extreme'};
peril_IDs       = {'FL_depth_monsoon' 'FL_duration_monsoon' 'FL_depth_cyclone' 'FL_duration_cyclone'};

% init
measures_impact = [];

% load flood entity, is always the same and does not need to be reloaded in
% the loop
entity = barisal_get_entity(entity_ref_year,'flood',entity_files);
entity.measures.cost = (entity.measures.cost)*0 + 10^4;

% loop over the four different flood perils
for peril_i = 1:length(peril_IDs)
    
    % load hazard
    hazard = barisal_get_hazard(hazard_ref_year,cc_scen,peril_IDs(peril_i),hazard_files);
    
    % calculate measures impact
    measures_impact_peril = climada_measures_impact_advanced(entity,hazard,'no');
    
    % combine measures_impact with measures_impact, that will contain the
    % added up benefits and cost_benefit_ratio
    if isempty(measures_impact)
        fprintf('\t-Start with measures impact from %s. \n', peril_IDs{peril_i})
        measures_impact = measures_impact_peril;
    else
        fprintf('\t-Add measures impact from %s. \n', peril_IDs{peril_i})
        measures_impact = climada_measures_impact_combine(measures_impact,measures_impact_peril);
    end
end
% finally create figure (multi peril adaptation cost curve)
figure
barisal_adaptation_cost_curve(measures_impact,[],[],[],0,0,1,0)
titlestr = sprintf('Barisal, %d, %s climate change', hazard_ref_year, CC_SCEN);
title({titlestr;'All perils combined (FL depth, duration, monsoon and cyclone, except TC wind)'})
hazard_name = 'all_perils';
entity_name = '2030_moderate_cc';
print(gcf,'-dpng',[results_dir filesep 'CBA_Barisal_BCC_' hazard_name '_' entity_name '.png']) % save as png


%% measures impact for specific scenario

% calculate the impact of measures for 1 scenario specified below

hazard_ref_year = 2014;
entity_ref_year = 2030;
peril_ID        = 'FL_depth_monsoon'; % Should be found in hazard filename
peril           = 'floods';           % Should be found in entity filename
cc_scen         = '';                 % Should be found in hazard filename

% retrieve relevant hazard and entity. The index (h_i/e_i) is for the cell
% arrays containing respective filenames
[hazard,h_i] = barisal_get_hazard(hazard_ref_year,cc_scen,peril_ID,hazard_files);
[entity,e_i] = barisal_get_entity(entity_ref_year,peril,entity_files);

% get name of entity/hazard
[~,hazard_name] = fileparts(hazard_files{h_i});
[~,entity_name] = fileparts(entity_files{e_i});

% calculate impact
fprintf('%s | %s | %d\n',entity.assets.comment,hazard.comment,hazard.reference_year) % to monitor
entity.measures.cost = (entity.measures.cost)*0 + 10^4;
measures_impact = climada_measures_impact_advanced(entity,hazard,'no');

% plot cost curve (barisal_adaptation_cost_curve is identical to climada_adaptation_cost_curve
% but sets currency to BDT)
% measures_impact.measures.cost = zeros(size(measures_impact.measures.cost))+10^4;
figure
barisal_adaptation_cost_curve(measures_impact,[],[],[],0,0,1,0)
print(gcf,'-dpng',[results_dir filesep 'CBA_' hazard_name '_' entity_name '.png']) % save as png

clear hazard_ref_year entity_ref_year peril_ID peril cc_scen hazard_name entity_name


%% stats for Gerbrand v Bork
ward_ndx = [32 34 27 28 33 29 30]; %shape index for wards 1-7
for ward_i = ward_ndx
    [POI.lon(find(ward_ndx == ward_i)),POI.lat(find(ward_ndx == ward_i))] ...
        = utm2ll_shift(mean(BCC_wards(ward_i).BoundingBox(:,1)),...
        mean(BCC_wards(ward_i).BoundingBox(:,2)));
    POI.name{find(ward_ndx == ward_i)}= BCC_wards(ward_i).UNION_NAME;
end
hazard = barisal_get_hazard(2030,'extreme','FL_duration_monsoon',hazard_files);
IFC = climada_hazard2IFC(hazard,POI,1);
clear POI ward_i ward_ndx
