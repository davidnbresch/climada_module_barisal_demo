%% Barisal Risk Calculations
clc
climada_global.waitbar = 0;
climada_global.EDS_at_centroid = 0;

%% Directories
barisal_data_dir= [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
entities_dir    = [barisal_data_dir filesep 'entities'];
hazards_dir     = [barisal_data_dir filesep 'hazards'];
results_dir     = [barisal_data_dir filesep 'results'];

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

%% entity

% This section constructs a cell array with entity file names for easy
% retrieval later with barisal_get_entity

% Entity excel file from Ecorys (unstructured)
entity_file_xls = [entities_dir filesep 'BCC_entity_170615_se_1.xls'];
% entity_file_xls = [entities_dir filesep 'Measure_Spatial_Zoning_se_1.xls']
% entity_file_xls = [entities_dir filesep 'Measure_Early_warning_system_se_1.xls'];
% entity_file_xls = [entities_dir filesep 'Meaure_Flood_proof_road_infrastructure_se_1.xls'];
% entity_file_xls = [entities_dir filesep 'BCC_entity_210615_se_1_population.xls'];

% Damage function file from Ecorys
% damfun_file_xls = [entities_dir filesep 'BCC_dmg_functions_170615.xls'];
damfun_file_xls = [entities_dir filesep 'BCC_dmg_functions_250615.xls'];
% Entity template file from global data dir
entity_temp_xls = [climada_global.data_dir filesep 'entities' filesep 'entity_template.xls'];
% Sheet names in Ecorys entity xls file
sheets          = {'Floods_2014' 'Floods_2030' 'Floods_2050' 'Cyclones_2014' 'Cyclones_2030' 'Cyclones_2050'};
% sheets          = {'2014' '2030' '2050'}; % for population entity

% Whether to re-read entity/damagefunctions from Ecorys xls (1) or load matfiles (0)
force_assets_re_read   = 0;
force_damfun_re_read   = 1;

for s_i = 1:length(sheets)
    clear entity
    
    % mat file name - separate mat file for each sheet in Ecorys entity xls
    entity_file_mat     = strrep(entity_file_xls,'.xls',['_' lower(sheets{s_i}) '.mat']);
    
    if exist(entity_file_mat,'file') && ~force_assets_re_read
        % Load mat file
        [~,fN] = fileparts(entity_file_mat);
        fprintf('entity %s already exists, skipping\n',fN)

        load(entity_file_mat)
        entity_files{s_i} = entity_file_mat;
        % append filename for consistency
        entity.assets.filename = entity_file_mat;
        % save corrected entity
        save(entity.assets.filename,'entity')  

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
        entity.assets = climada_xlsread('no',entity_file_xls,sheets{s_i},1);
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
        
        % set negative asset values to 0
        fprintf('setting deductible and cover\n')
        entity.assets.Deductible    = entity.assets.Value .* 0;
        entity.assets.Cover         = entity.assets.Value;
        
        % coord transformation from UTM to lat lon
        [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.lon, entity.assets.lat);
        
        % get reference year, comment from sheet name
        [~,entity.assets.reference_year]= str2num(strtok(sheets{s_i},'_'));
        entity.assets.comment           = strrep(sheets{s_i},'_',' ');
        entity.assets.filename          = entity_file_mat;
        
        fprintf('saving entity as %s\n',entity.assets.filename)
        save(entity.assets.filename,'entity')
        entity_files{s_i} = entity.assets.filename; % for filename consistency
    end
    
    if force_damfun_re_read
        load(entity_file_mat)
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
        save(entity.assets.filename,'entity')
        entity_files{s_i} = entity.assets.filename;
    end
end
clear entity_file_xls sheets ul_loc s_i tc_ndx damfun_file_xls entity_temp_xls 
clear fld_i flds nan_ndx entity_file_mat force_damfun_re_read force_assets_re_read

%% population entity & casualties

% Entity excel file from Ecorys (unstructured)
entity_file_xls = [entities_dir filesep 'BCC_entity_210615_se_1_population.xls'];

% Entity template file from global data dir
entity_temp_xls = [climada_global.data_dir filesep 'entities' filesep 'entity_template.xls'];
sheets          = {'2014' '2030' '2050'}; % for population entity

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
        entity = climada_entity_read_wo_assets(entity_temp_xls);
        entity.assets = climada_xlsread('no',entity_file_xls,sheets{s_i},1);
        entity.assets.comment = ['Population entity ' sheets{s_i}];
        entity.assets.reference_year = str2num(sheets{s_i});
        % coord transformation from UTM to lat lon
        [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.lon, entity.assets.lat);
        entity.assets.filename          = entity_file_mat; % for filename consistency
        fprintf('saving entity as %s\n',entity.assets.filename)
        save(entity.assets.filename,'entity')
        entity_files{s_i} = entity.assets.filename; 
    end
    
end
clear entity_file_xls sheets s_i entity_temp_xls entity_file_mat force_assets_re_read

year_i = 2014; climada_global.present_reference_year = year_i;
year_f = 2030; climada_global.future_reference_year  = year_f;
cc_scen = 'moderate'; ed_i = 0; clear EDS1 EDS2 EDS3 % init

for peril_ID = {'FL_depth_monsoon' 'FL_depth_cyclone'}
    ed_i = ed_i +1;
    % EDS1 for scenario hazard and entity in present reference year
    [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
    [entity,e_i] = barisal_get_entity(year_i,'population',entity_files);
    fprintf('***** EDS1 for %s | %s *****\n',char(entity.assets.comment),char(strtok(hazard.comment,',')))
    EDS1(ed_i)    = climada_EDS_calc(entity,hazard,'',1);
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
print(gcf,'-dpng',[results_dir filesep 'BCC_Barisal_population_waterfall_2030.png'])


%% measures construction

% structure measures, assign measures.hazard_event_set by searching through
% barisal entities folded and finding measure folders starting with the
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

% assign measures struct to entities and save each one again
measures                              = climada_measures_construct(measures,4);
measures.name{end-3}                  = 'Flood resilient crops';
measures.damagefunctions_map{end-3}   = '1320to1321';
measures.peril_ID{end-3}              = 'FL_';

% save measures struct (before adding entity file)
measures.filename = [entities_dir filesep 'BCC_measures_' datestr(now,'ddmmyy') '.mat'];
fprintf('saving measures as %s\n',measures.filename)
save(measures.filename,'measures')

for e_i = 1:length(entity_files)
    load(entity_files{e_i})
    
    % remove old measures, just in case
    if isfield(entity,'measures'),entity = rmfield(entity,'measures'); end
    
    measures.name{end-2}        = 'Spatial zoning';
    measures.entity_file{end-2} = [fileparts(entity_files{e_i}) filesep...
        'Measure_Spatial_Zoning_se_1_' lower(strrep(entity.assets.comment,' ','_')) '.mat'];
    
    measures.name{end-1}        = 'Early warning system';
    measures.entity_file{end-1} = [fileparts(entity_files{e_i}) filesep...
        'Measure_Early_warning_system_se_1_' lower(strrep(entity.assets.comment,' ','_')) '.mat'];
    
    measures.name{end}          = 'Flood proof road infrastructure';
    measures.entity_file{end}   = [fileparts(entity_files{e_i}) filesep...
        'Meaure_Flood_proof_road_infrastructure_se_1_' lower(strrep(entity.assets.comment,' ','_')) '.mat'];
    
    % drainage for correct refence year
    for m_i = 1:length(measures)
        if isempty(strfind(entity.assets.reference_year
        end
    end
    
    entity.measures = measures;
    save(entity.assets.filename,'entity')
end
% display measures in their respective colours
for measure_i = 1:length(measures.name),cprintf(measures.color_RGB(measure_i,:),'%s\n',measures.name{measure_i}); end
clear e_i e_dir_i e_dir_ e_subdir_ e_subdir_ m_i R G B m_i measure_i


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

cc_scen     = 'moderate'; %{'moderate' 'extreme'};
year_i      = 2014;
year_f      = [2030 2050];

EDS_only    = 1;

%init
% clear EDS1 EDS2 EDS3 EDS4 EDS5 measures_impact1 measures_impact3 measures_impact5

% for year_f = [2050 2030]
ed_i = 4;
for peril_ID = peril_IDs %peril_ID = peril_IDs(1)
    if strcmp(peril_ID,'TC')
        peril = 'cyclone';
    else
        peril = 'flood';
    end
    ed_i = ed_i+1; % counter for EDS entries
    
    year_f = 2030;
    climada_global.future_reference_year = year_f;
    
    % EDS1 for scenario hazard and entity in present reference year
    [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
    [entity,e_i] = barisal_get_entity(year_i,peril,entity_files);
    %EDS1    = barisal_get_EDS(EDS,entity_files{e_i},hazard_files{h_i});
    fprintf('\n***** EDS1 for %s | %s *****\n',char(entity.assets.comment),char(strtok(hazard.comment,',')))
    EDS1(ed_i)    = climada_EDS_calc(entity,hazard,'',1);
    scen_name1 = ['Today''s; expected damage; ' num2str(year_i)];
    fprintf('Annual expected damage: %2.2f mn\n',EDS1(ed_i).ED/1000000)
    if ~EDS_only
        measures_impact1(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
        measures_impact1(ed_i).filename = [results_dir filesep 'BCC_measures_impact_' num2str(year_i) '.mat'];
        for i = 1:length(measures_impact1(ed_i).EDS)
            % convert back to UTM
            [measures_impact1(ed_i).EDS(i).assets.X,measures_impact1(ed_i).EDS(i).assets.Y] = ...
                 ll2utm_shift(measures_impact1(ed_i).EDS(i).assets.lat,measures_impact1(ed_i).EDS(i).assets.lon);
        end
        save(measures_impact1(ed_i).filename,'measures_impact1')
    end
    
    % EDS2 for socio-economic growth scenario: present hazard, future entity 2030
    [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
    [entity,e_i] = barisal_get_entity(year_f,peril,entity_files);
    %EDS2    = barisal_get_EDS(EDS,entity_files{e_i},hazard_files{h_i});
    fprintf('\n***** EDS2 for %s | %s *****\n',char(entity.assets.comment),char(strtok(hazard.comment,',')))
    EDS2(ed_i)    = climada_EDS_calc(entity,hazard,'',1);
    scen_name2 = ['Increase; from economic; growth ' num2str(year_f)];
        fprintf('Annual expected damage: %2.2f mn\n',EDS2(ed_i).ED/1000000)
    
    % EDS3 for climate change scenario: future entity 2030, future hazard 2030
    [hazard,h_i] = barisal_get_hazard(year_f,cc_scen,peril_ID,hazard_files);
    [entity,e_i] = barisal_get_entity(year_f,peril,entity_files);
    fprintf('\n***** EDS3 for %s | %s *****\n',char(entity.assets.comment),char(strtok(hazard.comment,',')))
    EDS3(ed_i)    = climada_EDS_calc(entity,hazard,'',1);
    scen_name3 = ['Increase; from ' cc_scen '; climate change; ' num2str(year_f)];
    fprintf('Annual expected damage: %2.2f mn\n',EDS3(ed_i).ED/1000000)    
    if ~EDS_only
        measures_impact3(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
        measures_impact3(ed_i).filename = [results_dir filesep 'BCC_measures_impact_cc' cc_scen '_' num2str(year_f) '.mat'];
        for i = 1:length(measures_impact3(ed_i).EDS)
            % convert back to UTM
            [measures_impact3(ed_i).EDS(i).assets.X,measures_impact3(ed_i).EDS(i).assets.Y] = ...
                 ll2utm_shift(measures_impact3(ed_i).EDS(i).assets.lat,measures_impact3(ed_i).EDS(i).assets.lon);
        end
        save(measures_impact3(ed_i).filename,'measures_impact3')
    end
    
    year_f = 2050;
    climada_global.future_reference_year = year_f;
    
    % EDS4 for socio-economic growth scenario: present hazard, future entity 2050
    [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
    [entity,e_i] = barisal_get_entity(year_f,peril,entity_files);
    fprintf('\n***** EDS4 for %s | %s *****\n',char(entity.assets.comment),char(strtok(hazard.comment,',')))
    EDS4(ed_i)    = climada_EDS_calc(entity,hazard,'',1);
    scen_name4 = ['Increase; from economic; growth ' num2str(year_f)];
    fprintf('Annual expected damage: %2.2f mn\n',EDS4(ed_i).ED/1000000)
    
    % EDS5 for climate change scenario: future entity 2050, future hazard 2050
    [hazard,h_i] = barisal_get_hazard(year_f,cc_scen,peril_ID,hazard_files);
    [entity,e_i] = barisal_get_entity(year_f,peril,entity_files);
    fprintf('\n***** EDS5 for %s | %s *****\n',char(entity.assets.comment),char(strtok(hazard.comment,',')))
    EDS5(ed_i)    = climada_EDS_calc(entity,hazard,'',1);
    scen_name5 = ['Increase; from ' cc_scen '; climate change; ' num2str(year_f)];
    fprintf('Annual expected damage: %2.2f mn\n',EDS5(ed_i).ED/1000000)    
    if ~EDS_only
        measures_impact5(ed_i) = climada_measures_impact_advanced(entity,hazard,'no');
        measures_impact5(ed_i).filename = [results_dir filesep 'BCC_measures_impact_cc' cc_scen '_' num2str(year_f) '.mat'];
        for i = 1:length(measures_impact5(ed_i).EDS)
            % convert back to UTM
            [measures_impact5(ed_i).EDS(i).assets.X,measures_impact5(ed_i).EDS(i).assets.Y] = ...
                 ll2utm_shift(measures_impact5(ed_i).EDS(i).assets.lat,measures_impact5(ed_i).EDS(i).assets.lon);
        end
        save(measures_impact5(ed_i).filename,'measures_impact5')
    end
    
end  

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

% multi peril waterfall
fig = climada_waterfall_graph_multi_peril(0,'BDT',EDS1,scen_name1,EDS2,scen_name2,EDS3,scen_name3);
print(fig,'-dpng',[results_dir filesep 'BCC_waterfall_multi_peril_2030_' char(cc_scen) '.png'])
close

fig = climada_waterfall_graph_multi_peril(0,'BDT',EDS1,scen_name1,EDS4,scen_name4,EDS5,scen_name5);
print(fig,'-dpng',[results_dir filesep 'BCC_waterfall_multi_peril_2050_' char(cc_scen) '.png'])
close

%% measures impact report
EDS_report_xls = [results_dir filesep 'BCC_ED_report_measures_' datestr(now,'ddmmyy') '_test.xls'];

for m = [1 3 5]
    measures_impact = eval(['measures_impact' num2str(m)]);
    sheet = num2str(measures_impact(1).EDS(1).reference_year);
%     clear EDS_report fld
%     EDS_report = measures_impact(1).EDS(1:end-1); % control init
    sum_ED_at_centroid = 0;
    for k=1:length(EDS_report)
        EDS_report(k).ED_at_centroid    = zeros(size(EDS_report(k).ED_at_centroid)); % init
        fld(k).B_at_centroid            = zeros(size(EDS_report(k).ED_at_centroid)); % init
%         EDS_report(k).hazard            = measures_impact(1).hazard;
        EDS_report(k).peril_ID          = '';
    end
    for i = 1:length(measures_impact) % loop through each hazard/entity combo
        % ED at centroid w/o measure
        ctrl_ED_at_c = measures_impact(i).EDS(end).ED_at_centroid;

        missing = ~ismember({measures_impact(i).EDS(:).annotation_name},{EDS_report(:).annotation_name});
        if any(missing)
            fprintf('adding additional measure impacts to EDS_report:\n')
            for m = 1:length(measures_impact(i).EDS)
                fprintf('\t%s\n',measures_impact(i).EDS(m).annotation_name)
            end
            EDS_report = [EDS_report measures_impact(i).EDS(missing)];
        end
        
        for j = 1:length(EDS_report)
            ndx = find(strcmpi(EDS_report(j).annotation_name,{measures_impact(i).EDS(:).annotation_name}));
            if ~isempty(ndx)
                EDS_report(j).ED_at_centroid  = EDS_report(j).ED_at_centroid + measures_impact(i).EDS(ndx).ED_at_centroid;
                fld(j).B_at_centroid          = fld(j).B_at_centroid + (ctrl_ED_at_c -measures_impact(i).EDS(ndx).ED_at_centroid);

                EDS_report(ndx).hazard(end+1) = measures_impact(i).EDS(j).hazard;
                EDS_report(j).peril_ID      = [EDS_report(j).peril_ID measures_impact(i).EDS(ndx).peril_ID ' & '];
            else
                EDS_report(j).ED_at_centroid= EDS_report(j).ED_at_centroid+ctrl_ED_at_c;
                fld(j).B_at_centroid    = fld(j).B_at_centroid + 0;
            end
        end
    end
    for l = 1:length(EDS_report)
        EDS_report(l).B_at_centroid = fld(l).B_at_centroid;
    end
    % print to respective sheet in excel file
    climada_EDS_ED_at_centroid_report_xls(EDS_report,EDS_report_xls,sheet,'MI_at_centroid');
end
    

for ed_i = 1:length(EDS_report)
    climada_MB_plot(EDS_report(ed_i), 0,'BDT',0,0)
    shape_plotter(BCC_wards_ll,'','linewidth',1,'color',[81 81 81]/255);
    print(gcf,'-dpng',[results_dir filesep 'BCC_measure_benefit' strrep(EDS_report(ed_i).annotation_name,' ','_') '.png'])
    close
end  
      
% 
% % multi horizon waterfall
% for i = 1:5
%     eval(['EDS' num2str(i) '_combined = climada_EDS_combine(EDS' num2str(i) ')']);
% end
% fig = climada_waterfall_graph_2timehorizons(0,'BDT',EDS1,scen_name1,EDS4,scen_name4,EDS5,scen_name5);
% print(fig,'-dpng',[results_dir filesep 'BCC_waterfall_multi_peril_2050_' char(cc_scen) '.png'])
% close
% 
% % waterfall graph for today, eco, cc, 2030, eco 2050, cc 2050, tcr 2050
% fig = climada_waterfall_graph_2timehorizons([EDS1(1) EDS2(1) EDS3(1) EDS4(1) EDS5(1)],'AED');
% title(sprintf('Hazard type: %s','Monsoon'))
% foldername = sprintf('%sresults%sdamage_plots%sED_waterfall_from_%s.pdf', filesep,filesep,filesep,'Monsoon_2030_2050');
% print(fig,'-dpdf',[module_data_dir foldername])

% end
%% baseline EDS for ED report for Ecorys
% EDS_baseline = [EDS1 EDS2 EDS3 EDS4 EDS5]; save([results_dir filesep 'EDS_baseline.mat'],'EDS_baseline')
EDS_baseline = [EDS1 EDS2 EDS3]; %save([results_dir filesep 'EDS_baseline.mat'],'EDS_baseline')

% EDS_report_xls = [results_dir filesep 'BCC_ED_report_spatial_zoning_' datestr(now,'ddmmyy') '.xls'];
EDS_report_xls = [results_dir filesep 'BCC_ED_report_population_' datestr(now,'ddmmyy') '.xls'];

if exist(EDS_report_xls,'file'), delete(EDS_report_xls); end
% convert back to UTM
for ed_i = 1:length(EDS_baseline)
    [EDS_baseline(ed_i).assets.X,EDS_baseline(ed_i).assets.Y] = ...
        ll2utm_shift(EDS_baseline(ed_i).assets.lat,EDS_baseline(ed_i).assets.lon);
end

% write to excel
climada_EDS_ED_at_centroid_report_xls(EDS_baseline,barisal_get_entity(2014,'population',entity_files),EDS_report_xls);

clear peril_ID peril year_i year_f EDS_only EDS_report_xls h_i e_i peril_IDs CC_SCEN YEAR_F

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
