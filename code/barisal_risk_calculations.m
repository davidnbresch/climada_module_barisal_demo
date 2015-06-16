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
                if strcmp(hazard.peril_ID,'TC')
                    hazard.comment = strrep(strrep(hazard_file,'_',' '),'.mat','');
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
force_damfun_re_read   = 0;

for s_i = 1:length(sheets)
    clear entity
    
    entity_file_mat     = strrep(entity_file_xls,'.xls',['_' lower(sheets{s_i}) '.mat']);
    
    if exist(entity_file_mat,'file') && ~force_assets_re_read
        fprintf('entity %s .mat file already exists, skipping\n',lower(sheets{s_i}))
        entity_files{s_i} = entity_file_mat;
        load(entity_file_mat)
        entity.assets.filename = entity_file_mat;
        save(entity.assets.filename,'entity')
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

%% damage calc
EDS_save_file = [results_dir filesep 'BCC_EDS_' datestr(now,'ddmmyy') '.mat'];
EDS_load_file = [results_dir filesep 'BCC_EDS_090615.mat'];

EDS_force_recalc = 0;

if exist(EDS_load_file,'file') && ~EDS_force_recalc
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
clear max_val ndx EDS_force_recalc

%% damage & hazard plotter
% for ed_i = 1:length(EDS)
%     %     climada_ED_plot_per_point(EDS_(ed_i),BCC_wards);
%     climada_ED_plot(EDS(ed_i), 0,'BDT',0,0)
%     shape_plotter(BCC_wards_ll,'UNION_NAME','color',[81 81 81]/255,'linewidth',1.5)
%     fN = [results_dir filesep 'damage_plots' filesep ...
%         strrep(strrep(EDS(ed_i).annotation_name,' ','_'),'hazard','ED')];
%     print(gcf,'-dpng',[fN '.png'])
%     close
%
%     load(EDS(ed_i).hazard.filename)
%     figure
%     climada_hazard_plot_hr(hazard,0,[],[],0,0);
%     shape_plotter(BCC_wards_ll,'UNION_NAME','color',[81 81 81]/255,'linewidth',1.2)
%     print(gcf,'-dpng',[results_dir filesep 'hazard_plots' filesep ...
%         strrep(EDS(ed_i).annotation_name,' ','_') '.png'])
%     close
% end
%
% clear ed_i

%% measures construction
e_dir_ = dir(entities_dir);

%init
measures.asci_file = {}; measures.name={};measures.color={};measures.color_RGB = [];
measures.cost= [];measures.hazard_intensity_impact=[];
measures.hazard_high_frequency_cutoff=[];measures.hazard_event_set={};
measures.MDD_impact_a= [];measures.MDD_impact_b= [];measures.PAA_impact_a= [];
measures.PAA_impact_b= [];measures.damagefunctions_map={};
measures.risk_transfer_attachement = [];measures.risk_transfer_cover = [];
measures.peril_ID={}; measures.hazard_event_set_operator={};

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
                measures.hazard_event_set_operator{end+1} = 't';
            else
                measures.peril_ID{end+1} = 'FL';
                measures.hazard_event_set_operator{end+1} = 'p';
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
    
    if isfield(entity,'measures'),entity = rmfield(entity,'measures'); end
    entity.measures = measures;
    save(entity.assets.filename,'entity')
end

%%

% hazard_ref_year = 2030;
% entity_ref_year = 2030;
% peril_ID        = 'FL_depth_cyclone';
% peril           = 'floods';
% cc_scen         = 'moderate';
%
% [hazard,h_i] = barisal_get_hazard(hazard_ref_year,cc_scen,peril_ID,hazard_files);
% [entity,e_i] = barisal_get_entity(entity_ref_year,peril,entity_files);
%
% [~,hazard_name] = fileparts(hazard_files{h_i});
% [~,entity_name] = fileparts(entity_files{e_i});
%
% fprintf('%s | %s | %d\n',entity.assets.comment,hazard.comment,hazard.reference_year)
% measures_impact = climada_measures_impact_advanced(entity,hazard,'no');
% figure
% climada_adaptation_cost_curve(measures_impact)
% figure
% barisal_adaptation_cost_curve(measures_impact,[],[],[],0,0,0,0)
% print(gcf,'-dpng',[results_dir filesep 'CBA_' hazard_name '_' entity_name '.png'])
%
% clear hazard_ref_year entity_ref_year peril_ID peril cc_scen

%%
peril_IDs   = {'FL_depth_cyclone' 'FL_depth_monsoon' 'FL_duration_cyclone' 'FL_duration_monsoon','TC'};
cc_scen     = 'moderate'; %{'moderate' 'extreme'};
year_i      = 2014;
year_f      = 2030;

for cc_scen = {'moderate' 'extreme'}
    for year_f = [2030 2050]
        ed_i = 0;
        for peril_ID = peril_IDs
            if strcmp(peril_ID,'TC')
                peril = 'cyclone';
            else
                peril = 'flood';
            end
            ed_i = ed_i+1;
            [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
            [entity,e_i] = barisal_get_entity(year_i,peril,entity_files);
            %EDS1    = barisal_get_EDS(EDS,entity_files{e_i},hazard_files{h_i});
            fprintf('***** EDS1 for %s | %s %d *****\n',entity.assets.comment,strtok(hazard.comment,','),hazard.reference_year)
            EDS1(ed_i)    = climada_EDS_calc(entity,hazard);
            if ~strcmp(peril_ID,'TC')
                measures_impact = climada_measures_impact_advanced(entity,hazard,'no');
                measures_impact.filename = strrep(measures_impact.filename,'measures','measuresimpact');
                save(measures_impact.filename,'measures_impact')
                barisal_adaptation_cost_curve(measures_impact, [],[],[],[],0)
                plot_save_name =strrep(strrep(measures_impact.title_str,' ','_'),'|','-');
                print(gcf,'-dpng',[results_dir filesep plot_save_name  '.png'])
                close
            end
            
            [hazard,h_i] = barisal_get_hazard(year_i,'',peril_ID,hazard_files);
            [entity,e_i] = barisal_get_entity(year_f,peril,entity_files);
            %EDS2    = barisal_get_EDS(EDS,entity_files{e_i},hazard_files{h_i});
            fprintf('***** EDS2 for %s | %s *****\n',entity.assets.comment,hazard.comment)
            EDS2(ed_i)    = climada_EDS_calc(entity,hazard);
            if ~strcmp(peril_ID,'TC')
                measures_impact = climada_measures_impact_advanced(entity,hazard,'no');
                measures_impact.filename = strrep(measures_impact.filename,'measures','measuresimpact');
                save(measures_impact.filename,'measures_impact')
                barisal_adaptation_cost_curve(measures_impact, [],[],[],[],0)
                plot_save_name =strrep(strrep(measures_impact.title_str,' ','_'),'|','-');
                print(gcf,'-dpng',[results_dir filesep plot_save_name  '.png'])
                close
            end
            
            [hazard,h_i] = barisal_get_hazard(year_f,cc_scen,peril_ID,hazard_files);
            [entity,e_i] = barisal_get_entity(year_f,peril,entity_files);
            %EDS3    = barisal_get_EDS(EDS,entity_files{e_i},hazard_files{h_i});
            fprintf('***** EDS3 for %s | %s *****\n',entity.assets.comment,hazard.comment)
            EDS3(ed_i)    = climada_EDS_calc(entity,hazard);
            if ~strcmp(peril_ID,'TC')
                measures_impact = climada_measures_impact_advanced(entity,hazard,'no');
                measures_impact.filename = strrep(measures_impact.filename,'measures','measuresimpact');
                save(measures_impact.filename,'measures_impact')
                barisal_adaptation_cost_curve(measures_impact, [],[],[],[],0)
                plot_save_name =strrep(strrep(measures_impact.title_str,' ','_'),'|','-');
                print(gcf,'-dpng',[results_dir filesep plot_save_name  '.png'])
                close
            end            %                 climada_waterfall_graph(EDS1,EDS2,EDS3,'AED');
            %                 print(gcf,'-dpng',[results_dir filesep 'BCC_CBA_' ...
            %                     peril_ID '_' num2str(year_f) '_' cc_scen '.png'])
            %                 close
        end
        
        climada_waterfall_graph_multi_peril(0,'BDT',EDS1,EDS2,EDS3)
        print(gcf,'-dpng',[results_dir filesep 'BCC_CBA_multi_peril_' num2str(year_f) '_' cc_scen '.png'])
        close
    end
end

clear peril_ID peril year_i year_f
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
