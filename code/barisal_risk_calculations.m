%% Barisal Risk Calculations
climada_global.waitbar = 0;

%% Directories
barisal_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
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
hazard_file_tmp = 'Barisal_BCC_hazard_PIDSPEC_CCSCEN.mat';

CCSCEN  = {'2014' 'cc_2030_moderate' 'cc_2030_extreme' 'cc_2050_moderate' 'cc_2050_extreme'};
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
                if ~isempty(strfind(hazard_file,'duration'))
                    hazard.peril_ID = 'FL_';
                end
                hazard.filename = hazard_files{file_i};
                save(hazard.filename,'hazard')
            end
        end
    end
end
clear cc CCSCEN pid PID spec SPEC file_i hazard_file_tmp hazard_file

%% entity
entity_file_xls = [entities_dir filesep 'BCC_entity_030615_se_1.xls'];
damfun_file_xls = [entities_dir filesep 'BCC_dmg_functions_030615.xls'];
entity_temp_xls = [climada_global.data_dir filesep 'entities' filesep 'entity_template.xls'];


sheets          = {'Floods_2014' 'Floods_2030' 'Floods_2050' 'Cyclones_2014' 'Cyclones_2030' 'Cyclones_2050'};

force_assets_re_read   = 0;
force_damfun_re_read   = 1;

for s_i = 1:length(sheets)
    clear entity
    
    entity_file_mat     = strrep(entity_file_xls,'.xls',['_' lower(sheets{s_i}) '.mat']);
    
    if exist(entity_file_mat,'file') && ~force_assets_re_read
        fprintf('entity %s .mat file already exists, skipping\n',lower(sheets{s_i}))
        entity_files{s_i} = entity_file_mat;
        load(entity_file_mat)
        %continue
    else
        
        % assets
        entity.assets = climada_xlsread('no',entity_file_xls,sheets{s_i},1);
        entity        = barisal_entity_pre_process(entity);
        nan_ndx = isnan(entity.assets.lon);
        entity.assets.lon   (nan_ndx) = [];
        entity.assets.lat   (nan_ndx) = [];
        entity.assets.Value (nan_ndx) = [];
        entity.assets.Value(isnan(entity.assets.Value)) = 0;
        
        [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.lon, entity.assets.lat);
        ul_loc                          = strfind(sheets{s_i},'_');
        entity.assets.reference_year    = str2double(sheets{s_i}(ul_loc+1:end));
        entity.assets.comment           = strrep(sheets{s_i},'_',' ');
        entity.assets.filename          = entity_file_mat;
        
        fprintf('saving entity as %s\n',entity.assets.filename)
        save(entity.assets.filename,'entity')
        entity_files{s_i} = entity.assets.filename;
    end
    
    if force_damfun_re_read
        % damagefunctions
        entity.damagefunctions          = climada_xlsread(0,damfun_file_xls,'formatted',1);
        tc_ndx                          = strcmp(entity.damagefunctions.peril_ID,'TC');
        % widnspeed conversion from km/h to m/s
        entity.damagefunctions.Intensity(tc_ndx) = entity.damagefunctions.Intensity(tc_ndx)./3.6;
        entity.damagefunctions.units    (tc_ndx) = {'m/s'};
        
        entity.damagefunctions.MDD = sqrt(entity.damagefunctions.MDR);
        entity.damagefunctions.PAA = sqrt(entity.damagefunctions.MDR);
        
        % discount
        entity.discount          = climada_xlsread(0,entity_temp_xls,'discount',1);
        
        fprintf('saving entity as %s\n',entity.assets.filename)
        save(entity.assets.filename,'entity')
        entity_files{s_i} = entity.assets.filename;
    end
end
clear entity_file_xls sheets ul_loc s_i tc_ndx damfun_file_xls

%% entity filename
% entity_file_tmp = 'Spreadsheet 100x100 Assets at risk PID1 060515PID2.mat';
%
% PID1 = {'Cyclones' 'Flooding'};
% PID2 = {'_cyclone_windspeed' '_flood_duration' '_flood_depth'};
%
% file_i = 0; entity_files = {};
% for pid1 = PID1
%     for pid2 = PID2
%         entity_file = entity_file_tmp;
%         entity_file = strrep(entity_file,'PID1',char(pid1));
%         entity_file = strrep(entity_file,'PID2',char(pid2));
%         if exist([entities_dir filesep entity_file],'file')
%             file_i = file_i+1;
%             entity_files{file_i} = [entities_dir filesep entity_file];
%         end
%     end
% end
% clear pid1 PID1 pid2 PID2 file_i entity_file_tmp entity_file


%% damage calc
EDS_save_file = [results_dir filesep 'BCC_EDS_' datestr(now,'ddmmyy') '.mat'];
EDS_load_file = [results_dir filesep 'BCC_EDS_090615.mat'];
if exist(EDS_load_file,'file')
    load(EDS_load_file)
else
    EDS        = climada_EDS_multi_calc(entity_files,hazard_files,EDS_save_file,1,0);
end

for ed_i = 1:length(EDS)
    ndx = find([EDS.reference_year] == EDS(ed_i).reference_year);
    max_val = [];
    for ndx_i = ndx
        max_val = max([max_val EDS(ndx_i).assets.Value],[],2);
    end
    
    EDS(ed_i).Value_total = sum(max_val);
end
clear max_val ndx
%% plotter
for ed_i = 1:length(EDS)
    % %     climada_ED_plot_per_point(EDS_(ed_i),BCC_wards);
    %     climada_ED_plot(EDS(ed_i), 0,'BDT',0,0)
    %     shape_plotter(BCC_wards_ll,'UNION_NAME','color',[81 81 81]/255,'linewidth',1.5)
    %     fN = [results_dir filesep 'damage_plots' filesep ...
    %         strrep(strrep(EDS(ed_i).annotation_name,' ','_'),'hazard','ED')];
    %     print(gcf,'-dpng',[fN '.png'])
    % %     save([fN '.asc'],,'-ascii','-double'
    %     close
    
    load(EDS(ed_i).hazard.filename)
    figure
    climada_hazard_plot_hr(hazard,0,[],[],0,0);
    shape_plotter(BCC_wards_ll,'UNION_NAME','color',[81 81 81]/255,'linewidth',1.2)
    print(gcf,'-dpng',[results_dir filesep 'hazard_plots' filesep ...
        strrep(EDS(ed_i).annotation_name,' ','_') '.png'])
    close
end

clear ed_i

%% hazard - entity pairs
%
% hazard_entity   = {};
% flood_type      = {'duration' 'depth'};
%
% for h_i =1:length(hazard_files)
%     [h_fP,h_fN,h_fE] = fileparts(hazard_files{h_i});
%     load(hazard_files{h_i})
%     for e_i =1:length(entity_files)
%         [e_fP,e_fN,e_fE] = fileparts(entity_files{e_i});
%         load(entity_files{e_i})
%         if hazard.reference_year==entity.assets.reference_year
%             && ~isempty(strfind(h_fN,'FL'))
%             && ~isempty(strfind(e_fN,'floods'))
%             for ft_i = 1:length(flood_type)
%                 if ~isempty(strfind(hazard_files{h_i},flood_type{ft_i}))
%
%         end
%     end
% end





%% measures impact
e_dir_ = dir(entities_dir);

%init
measures.asci_file = {}; measures.name={};measures.color={};measures.color_RGB = [];
measures.cost= [];measures.hazard_intensity_impact=[];
measures.hazard_high_frequency_cutoff=[];measures.hazard_event_set={};
measures.MDD_impact_a= [];measures.MDD_impact_b= [];measures.PAA_impact_a= [];
measures.PAA_impact_b= [];measures.damagefunctions_map={};
measures.risk_transfer_attachement = [];measures.risk_transfer_cover = [];
measures.peril_ID={};

for e_dir_i = 1:length(e_dir_)
    if e_dir_(e_dir_i).isdir && ...
            ~isempty(strfind(upper(e_dir_(e_dir_i).name),'MEASURES'))
        
        e_subdir_ = dir([entities_dir filesep e_dir_(e_dir_i).name]);
        
        for m_i = 1:length(e_subdir_)
            if e_subdir_(m_i).name(1) == '.', continue; end
            measures.asci_file{end+1} ...
                = [entities_dir filesep e_dir_(e_dir_i).name filesep e_subdir_(m_i).name];
            name = strrep([e_dir_(e_dir_i).name],'Measures_','');
            name(1) = upper(name(1)); name = strrep(name,'_',' ');
            
            measures.name{end+1} = name;
            R = rand; G = rand; B = rand;
            measures.color{end+1} = [num2str(R) ' ' num2str(G) ' ' num2str(B)];
            measures.color_RGB(end+1,:) = [R; G; B];
            measures.cost(end+1)                            = [1];
            measures.hazard_intensity_impact(end+1)         = [0];
            measures.hazard_high_frequency_cutoff(end+1)    = [0];
            measures.hazard_event_set{end+1}                =  measures.asci_file{end};
            measures.MDD_impact_a(end+1)                    = [1];
            measures.MDD_impact_b(end+1)                    = [0];
            measures.PAA_impact_a(end+1)                    = [1];
            measures.PAA_impact_b(end+1)                    = [0];
            measures.damagefunctions_map{end+1}             = 'nil';
            measures.risk_transfer_attachement(end+1)       = [0];
            measures.risk_transfer_cover(end+1)             = [0];
            if ~isempty(strfind(e_subdir_(m_i).name,'Duration'))
                measures.peril_ID{end+1} = 'FL_';
            else
                measures.peril_ID{end+1} = 'FL';
            end
            fprintf('%s \t%s\n',name,measures.peril_ID{end})
        end
    end
end

measures.filename = [entities_dir filesep 'BCC_measures_' datestr(now,'ddmmyy') '.mat'];
fprintf('saving measures as %s\n',measures.filename)
save(measures.filename,'measures')

for e_i = 1:length(entity_files)
    if isempty(strfind(entity_files{1},'floods')),continue;end
    load(entity_files{e_i})
    
    % discount
    entity.discount          = climada_xlsread(0,entity_temp_xls,'discount',1);
    
    if isfield(entity,'measures'),entity = rmfield(entity,'measures'); end
    entity.measures = measures;
    save(entity.assets.filename,'entity')
end

% for h_i = 1:length(hazard_files)
%     if isempty(strfind(entity_files{1},'FL')),continue;end
%     load(hazard_files{h_i})
%
%     for e_i = 1:length(entity_files)
%         if isempty(strfind(entity_files{1},'floods')),continue;end
%         load(entity_files{e_i})
%         if isfield(entity,'measures'),entity = rmfield(entity,'measures'); end
%         entity.measures = measures;
%         for m_i = 1:length(measures.name)
%             if isempty(strfind(entity_files{1},'FL')),continue;end
%
%             hazard_mod_file = strrep(hazard_files{h_i},'.mat',['_' strrep(lower(measures.name),'.asc','') '.mat']);
%             climada_distributed_measures(measures(m_i).asci_file,hazard,hazard_mod_file)
%             measures.hazard_event_set{m_i} = hazard_mod_file;
%         end
%     end
% end
%
% hazard_w_measures = climada_distributed_measures(







