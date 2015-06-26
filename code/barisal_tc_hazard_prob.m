
%% plot historical tc tracks
tc_track_file = [climada_global.data_dir filesep 'tc_tracks' filesep 'tracks_n_indian_proc'];
load(tc_track_file)

fig = climada_figuresize(0.5,0.7);
check_country = 'Bangladesh';
keep_boundary = 0;
climada_plot_world_borders(1,check_country,'',keep_boundary,'');
track_count = numel(tc_track);
for track_i = 1:track_count
    h = climada_plot_tc_track_stormcategory(tc_track(track_i),5,[]);
end
%add legend, makes it quite slow
climada_plot_tc_track_stormcategory(0,8,1);
axis equal
axis([79 98 14 28])
% axis([70 110 06 32])   
filename = [filesep 'results' filesep 'TC_tracks_NIO_hist.pdf'];
print(fig,'-dpdf',[climada_global.data_dir filename])
fprintf('figure saved in %s \n', filename) 
datestr(tc_track(1).datenum(1))
datestr(tc_track(end).datenum(end))



%% create tc track probabilistic for barisal
tc_track_file = [climada_global.data_dir filesep 'tc_tracks' filesep 'tracks_n_indian_proc'];
load(tc_track_file)

tc_track_save = [tc_track_file '_prob.mat'];
ens_size      = 9;
ens_amp       = 0.2; %degree
maxangle      = pi/4;
tc_track      = climada_tc_random_walk_position_windspeed(tc_track,tc_track_save,ens_size,ens_amp,maxangle,1, 0);
% ens_amp  = [];
% Maxangle = [];
% tc_track_out  = climada_tc_random_walk(tc_track,ens_size,ens_amp,Maxangle,0);
save(tc_track_save, 'tc_track')
load(tc_track_save)


%% create tc track figure
% load BCC boundaries file
shp_mat_file = [climada_global.data_dir filesep 'results' filesep 'BCC_boundary_shp.mat'];
load(shp_mat_file)
fig = climada_figuresize(0.5,0.7);
event_i = 173;
check_country = 'Bangladesh';
keep_boundary = 0;
climada_plot_world_borders(1,check_country,'',keep_boundary,'');
for t_i = (event_i-1)*(ens_size+1)+1:1:(event_i+0)*(ens_size+1)  %1:numel(tc_track)
    if tc_track(t_i).orig_event_flag == 1
        plot(tc_track(t_i).lon, tc_track(t_i).lat,'.-r','markersize',3)
        hold on
    else
        plot(tc_track(t_i).lon, tc_track(t_i).lat,'.-b','markersize',3)
        hold on
    end
end
for s_i = 1:numel(BCC_boundary)
    plot(BCC_boundary(s_i).X, BCC_boundary(s_i).Y,'-k');
end
axis([82 97 17 28])
axis equal
% axis([82 97 17 28])
axis([70 110 06 32])
xlabel('Longitude'); ylabel('Latitude')
foldername  = [filesep 'results' filesep 'Sidr_probabilistic_daughters.pdf'];
print(fig,'-dpdf',[climada_global.data_dir foldername])

% datestr(hazard.datenum(172*4+1:173*4))
% datestr(hazard.datenum((event_i-1)*(ens_size+1)+1:1:event_i*(ens_size+1)))

climada_plot_probabilistic_wind_speed_map(tc_track, (event_i-1)*(ens_size+1)+1)

%% footprint figure
% climada_plot_tc_footprint(hazard,tc_track((event_i-1)*(ens_size+1)+1))
% caxis_range = '';
% res=climada_hazard_plot(hazard,(event_i-1)*(ens_size+1)+1);
% 
% load([climada_global.modules_dir filesep 'barisal_demo' filesep 'data' filesep 'entities' filesep 'Barisal_BCC_1km_100.mat'])
% focus_region = [70 110 06 32];
% check_mode = '';
% tc_track_1 = tc_track((event_i-1)*(ens_size+1)+1);
% hazard = climada_event_damage_data_tc(tc_track_1,entity,'',0,check_mode,focus_region);
% climada_event_damage_animation


%% create probabilistic tc hazard set
centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
load(centroids_file)
hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_prob'];
hazard = climada_tc_hazard_set(tc_track, hazard_set_file, centroids);

% tweek the frequencies
hazard.frequency_ori = hazard.frequency;
hazard.frequency     = hazard.frequency_ori*6;
% ori_flag          = logical(hazard.orig_event_flag);
% hazard.frequency(ori_flag) = hazard.frequency_ori(ori_flag);
save(hazard_set_file,'hazard')

% add tc track category manually
% hazard.category = [tc_track.category];
% hazard.category(hazard.category<0) = 0;

% just loading not calculating
hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_prob'];
load(hazard_set_file)

%% view wind results in Barisal (centroid ID 30)
centroid_ID = [30];
IFC = climada_hazard2IFC(hazard, centroid_ID);
close all
climada_IFC_plot(IFC,0)
foldername  = [filesep 'results' filesep 'TC_wind_intensity_Barisal.pdf'];
print(gcf,'-dpdf',[climada_global.data_dir foldername])


%% was the max event Sidr? (Nov 2007)
orig_event_flag = logical(hazard.orig_event_flag);
[int_max indx]  = max(full(hazard.intensity(orig_event_flag,centroid_ID)));
indx = (indx-1)*(ens_size+1)+1;
hazard.name{indx}
datestr(hazard.datenum(indx))
% yes, this is Sidr: 46.8 m/s at centroid_ID 30, 10 Nov 2007

% datestr(hazard.datenum(172*4+1:173*4))
% int_Sidr = full(hazard.intensity(172*10+1,centroid_ID))

% event_i = 172*4+1;
% figure
% res=climada_hazard_plot(hazard,event_i,'','','');
% res=climada_hazard_plot(hazard,event_i+1,'','','');
% figure
% res=climada_hazard_plot(hazard,event_i+2,'','','');
% figure
% res=climada_hazard_plot(hazard,event_i+3,'','','');



















