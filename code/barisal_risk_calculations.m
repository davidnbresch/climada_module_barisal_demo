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
BCC_wards_ll = BCC_wards;
[BCC_wards_ll.x] = BCC_wards_ll.X; BCC_wards_ll = rmfield(BCC_wards_ll,'X');
[BCC_wards_ll.y] = BCC_wards_ll.Y; BCC_wards_ll = rmfield(BCC_wards_ll,'Y');
[BCC_wards_ll.X] = BCC_wards_ll.lon; BCC_wards_ll = rmfield(BCC_wards_ll,'lon');
[BCC_wards_ll.Y]= BCC_wards_ll.lat; BCC_wards_ll = rmfield(BCC_wards_ll,'lat');

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

%% entity filename
entity_file_tmp = 'Spreadsheet 100x100 Assets at risk PID1 060515PID2.mat';

PID1 = {'Cyclones' 'Flooding'};
PID2 = {'_cyclone_windspeed' '_flood_duration' '_flood_depth'};

file_i = 0; entity_files = {};
for pid1 = PID1
    for pid2 = PID2
        entity_file = entity_file_tmp;
        entity_file = strrep(entity_file,'PID1',char(pid1));
        entity_file = strrep(entity_file,'PID2',char(pid2));
        if exist([entities_dir filesep entity_file],'file')
            file_i = file_i+1;
            entity_files{file_i} = [entities_dir filesep entity_file];
        end
    end
end
clear pid1 PID1 pid2 PID2 file_i entity_file_tmp entity_file

%% entity modifier
% fmt = '%s';
% for file_i = 1:length(entity_files)
%     [~, fN, ~] = fileparts(entity_files{file_i});
%     ld_msg  = sprintf('loading %s',fN);
%     fprintf(fmt, ld_msg);   fmt     = [repmat('\b',size(ld_msg)) '%s']; 
%     load(entity_files{file_i});
%     
%     if isfield(entity.assets,'Category')
%         ndx = strcmp(entity.assets.Category,'Agriculture_Crops');
%         entity.assets.Value     (ndx)   = entity.assets.Value       (ndx)./20;
%         entity.assets.Value_2030(ndx)   = entity.assets.Value_2030  (ndx)./20;
%         entity.assets.Value_2050(ndx)   = entity.assets.Value_2050  (ndx)./20;
%     
%         sv_msg  = sprintf('saving %s',fN);
%         fprintf(fmt,sv_msg);    fmt     = [repmat('\b',size(sv_msg)) '%s']; 
%         save(entity_files{file_i},'entity');
%     end
% end
% fprintf('\n'); clear ndx ld_msg ld_fmt sv_msg sv_fmt

%% damage calc
token = 'duration';

EDS_save_file = [results_dir filesep 'EDS_' token '.mat'];

EDS_ = climada_EDS_multi_calc(entity_files,hazard_files,EDS_save_file,1,token,0);

clear token
%% plotter
for ed_i = 1:length(EDS_)
%     climada_ED_plot_per_point(EDS_(ed_i),BCC_wnards);
    climada_ED_plot(EDS_(ed_i), 0,BCC_wards_ll,'UNION_NAME','BDT',6)
    print(gcf,'-dpdf',[results_dir filesep ...
        strrep(strrep(EDS_(ed_i).annotation_name,' ','_'),'hazard','ED') '.pdf'])
    pause(0.5)
    close
    
    load(EDS_(ed_i).hazard.filename)
    figure
    climada_hazard_plot_hr(hazard,0,[],[],0,0,BCC_wards_ll);
    print(gcf,'-dpdf',[results_dir filesep ...
        strrep(EDS_(ed_i).annotation_name,' ','_') '.pdf'])
    close
end

clear ed_i

%% difference (rain only)
cyclone_depth_ndx = [1 3 4 7 9 11];
climada_waterfall_graph(EDS_(1),EDS_(4),EDS_(3), 'AED');

monsoon_depth_ndx = [2 5 6 8 10 12];
climada_waterfall_graph(EDS_(2),EDS_(6),EDS_(5), 'AED');






