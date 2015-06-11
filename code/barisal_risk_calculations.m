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
SPEC= {'' '_rain_only'};

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
            end
        end
    end
end
clear cc CCSCEN pid PID spec SPEC file_i hazard_file_tmp hazard_file

%% entity
entity_file_xls = [entities_dir filesep 'BCC_entity_030615_se_1.xls'];
damfun_file_xls = [entities_dir filesep 'BCC_dmg_functions_030615.xls'];

sheets          = {'Floods_2014' 'Floods_2030' 'Floods_2050' 'Cyclones_2014' 'Cyclones_2030' 'Cyclones_2050'};

force_re_read   = 0;

for s_i = 1:length(sheets)
    clear entity
    
    entity_file_mat     = strrep(entity_file_xls,'.xls',['_' lower(sheets{s_i}) '.mat']);
    
    if exist(entity_file_mat,'file') && ~force_re_read
        fprintf('entity %s .mat file already exists, skipping\n',lower(sheets{s_i}))
        entity_files{s_i} = entity_file_mat;
        continue
    end
        
    % assets
    entity.assets = climada_xlsread('no',entity_file_xls,sheets{s_i},1);
    entity        = barisal_entity_pre_process(entity);
    [entity.assets.lon, entity.assets.lat] = utm2ll_shift(entity.assets.lon, entity.assets.lat);
    ul_loc                          = strfind(sheets{s_i},'_');
    entity.assets.reference_year    = str2double(sheets{s_i}(ul_loc+1:end));
    entity.assets.comment           = strrep(sheets{s_i},'_',' ');
    entity.assets.filename          = entity_file_mat;
    
    % damagefunctions
    entity.damagefunctions          = climada_xlsread(0,damfun_file_xls,'formatted',1);
    tc_ndx                          = strcmp(entity.damagefunctions.peril_ID,'TC');
    % widnspeed conversion from km/h to m/s
    entity.damagefunctions.Intensity(tc_ndx) = entity.damagefunctions.Intensity(tc_ndx)./3.6;
    entity.damagefunctions.units    (tc_ndx) = {'m/s'};
    
    fprintf('saving entity as %s\n',entity.assets.filename)
    save(entity.assets.filename,'entity')
    entity_files{s_i} = entity.assets.filename;
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

%% entity modifier
% fmt = '%s';
% for file_i = 1:length(entity_files)
%     [~, fN, ~] = fileparts(entity_files{file_i});
%     ld_msg  = sprintf('loading %s',fN);
%     fprintf(fmt, ld_msg);   fmt     = [repmat('\b',size(ld_msg)) '%s']; 
%     load(entity_files{file_i});
%     
%     % damagefunctions
%     entity.damagefunctions          = climada_xlsread(0,damfun_file_xls,'formatted',1);
%     tc_ndx                          = strcmp(entity.damagefunctions.peril_ID,'TC');
%     % widnspeed conversion from km/h to m/s
%     entity.damagefunctions.Intensity(tc_ndx) = entity.damagefunctions.Intensity(tc_ndx)./3.6;
%     entity.damagefunctions.units    (tc_ndx) = {'m/s'};
% 
%     entity.damagefunctions.MDD = sqrt(entity.damagefunctions.MDR);
%     entity.damagefunctions.PAA = sqrt(entity.damagefunctions.MDR);
% 
%     sv_msg  = sprintf('saving %s',fN);
%     fprintf(fmt,sv_msg);    fmt     = [repmat('\b',size(sv_msg)) '%s']; 
%     save(entity_files{file_i},'entity');
% 
%     
% %     if isfield(entity.assets,'Category')
% %         ndx = strcmp(entity.assets.Category,'Agriculture_Crops');
% %         entity.assets.Value     (ndx)   = entity.assets.Value       (ndx)./20;
% %         entity.assets.Value_2030(ndx)   = entity.assets.Value_2030  (ndx)./20;
% %         entity.assets.Value_2050(ndx)   = entity.assets.Value_2050  (ndx)./20;
% %     
% %         sv_msg  = sprintf('saving %s',fN);
% %         fprintf(fmt,sv_msg);    fmt     = [repmat('\b',size(sv_msg)) '%s']; 
% %         save(entity_files{file_i},'entity');
% %     end
% end
% fprintf('\n'); clear ndx ld_msg ld_fmt sv_msg sv_fmt

%% damage calc
EDS_save_file = [results_dir filesep 'BCC_EDS_' datestr(now,'ddmmyy') '.mat'];

if exist(EDS_save_file,'file')
    load(EDS_save_file)
else
    EDS        = climada_EDS_multi_calc(entity_files,hazard_files,EDS_save_file,1,0);
end
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







