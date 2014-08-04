function hazard=bangladesh_hazard_TEST(force_recalc)
% climada
% NAME:
%   bangladesh_hazard_TEST
% PURPOSE:
%   TEST the tropical cyclone (TC), torrential rain (TR) and storm surge 
%   (TS) raw hazard creation:
%   1) get centroids for the test country (eg Bangladesh, see PARAMETERS)
%      if they do not exist, try to run GDP_entity in order to create them
%   2) create TC wind hazard event set
%   3) create TR torrential rain hazard event set
%   4) create TS storm surge hazard event set
%   show the result
%
%   In essence, you define the country and the code checks the generation
%   of centroids, TC, TR and TS hazard event sets
%
%   in essence a caller for codes
%   - climada_tc_hazard_set
%   - climada_tr_hazard_set
%   - tc_surge_hazard_create
%
%   see e.g. tc_surge_plot_3d for 3D plots of surge fields
% CALLING SEQUENCE:
%   bangladesh_hazard_TEST(force_recalc)
% EXAMPLE:
%   bangladesh_hazard_TEST
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   force_recalc: if =1, recalculate the hazard sets, even if they exist
%       (good for TEST while editing the code)
% OUTPUTS:
%   writes a couple files, such as entity, assets, bathymetry and
%       hazard event sets. The hazard event set returned is the last one,
%       i.e. the TS storm surge set
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20140420
% David N. Bresch, david.bresch@gmail.com, 20140804, major revision
%-

hazard    = []; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('force_recalc','var'), force_recalc = 0;end

% PARAMETERS
%
use_cbfreeze=0; % default =1 after MATLAB version 8.1.0.604 (R2013a)
%
% in essence, only the TEST country, TEST location, TEST_probabilistic
% and the TC track set needs to be defined,
% all further parameters below should be working
TEST_country_name='Bangladesh'; tc_track_file='tracks.nio.txt';
%
TEST_location.name='  Barisal'; % first two spaces for nicer labeling
TEST_location.longitude=90+30/60+0/3600;
TEST_location.latitude =22+48/60+0/3600;
%TEST_location=''; % set TEST_location='' to omit labeling
%
% whether we run historic only (=0, default, good enough to test) or fully probabilistic
TEST_probabilistic = 0; % default=0, since fast to check
%
% see comment above, unlikely one needs to change parameters below
%
% 1) centroids for study region
% -----------------------------
% define the file with centroids (geo-locations of the points we later
% evaluate and store storm surge heights at)
% see climada_create_GDP_entity to create centroids file
% centroids_file = [climada_global.additional_dir filesep 'tc_surge' filesep ...
%     'data' filesep 'system'   filesep TEST_country_name '_centroids.mat'];
centroids_file = [climada_global.system_dir filesep 'centroids_' TEST_country_name '.mat'];
% if the centroids are generated in the present code, the entitity is also stored (not needed for this TEST)
entity_file    = [climada_global.additional_dir filesep 'tc_surge' filesep ...
    'data' filesep 'entities' filesep TEST_country_name '_assets.mat'];
%
% 2) tropical cyclone (TC) tracks
% -------------------------------
% set UNISYS TC track data file (for info, see climada_tc_read_unisys_database)
unisys_file    = [climada_global.additional_dir filesep 'barisal_demo' filesep ...
    'data' filesep 'tc_tracks' filesep tc_track_file];
%
% 3) bathymetry parameters set in tc_surge_hazard_create
%
% 4) surge hazard event set
% -------------------------
% define the hazard event set file to store the TEST hazard event set
hazard_set_file_tc = [climada_global.additional_dir filesep 'barisal_demo' filesep 'data' filesep 'hazards' filesep TEST_country_name '_hazard_TC.mat'];
hazard_set_file_ts = [climada_global.additional_dir filesep 'barisal_demo' filesep 'data' filesep 'hazards' filesep TEST_country_name '_hazard_TS.mat'];
hazard_set_file_tr = [climada_global.additional_dir filesep 'barisal_demo' filesep 'data' filesep 'hazards' filesep TEST_country_name '_hazard_TR.mat'];

% Calculations start
% ==================


% 1) read the centroids
% ---------------------
if exist(centroids_file,'file') %% && ~force_recalc
    load(centroids_file) % load centroids
else
    % invoke the GDP_entity moduke to generate centroids and entity
    TEST_country_name_tmp=TEST_country_name;
    if strcmp(TEST_country_name,'Vietnam')   , TEST_country_name_tmp='Viet Nam';end
    if strcmp(TEST_country_name,'ElSalvador'), TEST_country_name_tmp='El Salvador';end
    [centroids,entity]=climada_create_GDP_entity(TEST_country_name_tmp);
    save(centroids_file,'centroids');
    % note: assets are not needed for the TEST, but since it's convenient to store them for later use
    save(entity_file,'entity');
    % visualize assets on map
    climada_plot_entity_assets(entity,centroids,TEST_country_name);
end

% prep the region we need
centroids_rect = [min(centroids.Longitude) max(centroids.Longitude) min(centroids.Latitude) max(centroids.Latitude)];


% 2) HAZARD TROPICAL CYCLONE (TC) WIND: create TC hazard event set
% =================================================================
% note that if not in TEST mode, one might already have a fully
% probabilistic TC hazard event set hence does not need to (re)create
% in order for this TEST environment to work properly and almost
% independent of core climada, we (re)create the TC hazard event set here

if ~exist(hazard_set_file_tc,'file') || force_recalc
    
    % read tracks from unisys database file (txt)
    tc_track        = climada_tc_read_unisys_database(unisys_file);
    unisys_file_mat = strrep(strrep(unisys_file,'s.','s_'),'.txt','_hist.mat');
    save(unisys_file_mat, 'tc_track')
    
    if TEST_probabilistic
        
        if exist('climada_tc_track_wind_decay_calculate','file')
            % wind speed decay at track nodes after landfall
            [a p_rel]  = climada_tc_track_wind_decay_calculate(tc_track,1);
        else
            fprintf('NO inland decay, consider module tc_hazard_advanced\n');
        end
        
        tc_track = climada_tc_random_walk(tc_track); % overwrite
        
        if exist('climada_tc_track_wind_decay_calculate','file')
            % add the inland decay correction to all probabilistic nodes
            tc_track   = climada_tc_track_wind_decay(tc_track, p_rel,1);
        end
        
        % plot the tracks
        figure('Name','TC tracks','Color',[1 1 1]);
        hold on
        for event_i=1:length(tc_track) % plot all tracks
            plot(tc_track(event_i).lon,tc_track(event_i).lat,'-b');
        end % event_i
        % overlay historic (to make them visible, too)
        for event_i=1:length(tc_track)
            if tc_track(event_i).orig_event_flag
                plot(tc_track(event_i).lon,tc_track(event_i).lat,'-r');
            end
        end % event_i
        climada_plot_world_borders(2)
        box on
        axis equal
        axis(centroids_rect);
        xlabel('blue: probabilistic, red: historic');
        
        %save probabilistic track set
        save(strrep(unisys_file_mat,'_hist.mat','_prob.mat'),'tc_track')
        
    end
    
    % HAZARD TROPICAL WINDSTORMS: generate all the wind footprints and save hazard set
    hazard = climada_tc_hazard_set(tc_track, hazard_set_file_tc, centroids);
    
else
    fprintf('loading hazard TC wind from %s\n',hazard_set_file_tc);
    load(hazard_set_file_tc);
end

fprintf('TC: max(max(hazard.arr))=%f\n',full(max(max(hazard.arr)))); % a kind of easy check


% FIGURE: show biggest TC event
[~,max_tc_pos] = max(sum(hazard.arr,2)); % the maximum TC intensity

main_fig = climada_figuresize(0.75,0.8);
% main_fig = figure('Name','TEST','Position',[89 223 1014 413],'Color',[1 1 1]);
subplot(3,1,1)
values   = full(hazard.arr(max_tc_pos,:)); % get one TC footprint
centroids.Longitude   = hazard.lon; % as the gridding routine needs centroids
centroids.Latitude    = hazard.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(values,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
plot(centroids.Longitude,centroids.Latitude,'.r','MarkerSize',1);
if isfield(centroids,'onLand')
    water_points=find(centroids.onLand==0);
    plot(centroids.Longitude(water_points),centroids.Latitude(water_points),'.b','MarkerSize',1);
end
box on
climada_plot_world_borders
axis equal
axis(centroids_rect);
colorbar
title(sprintf('TC Wind Field [m/s] (Event %i)',max_tc_pos));
fprintf('max event %i\n\n',max_tc_pos);
cmap = climada_colormap(hazard.peril_ID);
colormap(cmap)
if use_cbfreeze
    freezeColors %freeze this plot's colormap
    cbfreeze(colorbar)
end
% up to here, hazard contains the tropical cyclone (TC) hazard event set


% 3) HAZARD TORRENTIAL RAIN
% ==========================
if ~exist(hazard_set_file_tr,'file') || force_recalc
    
    % placeholder for wind hazard
    hazard_TC = hazard;
    
    % load tc tracks if generated before
    if ~exist('tc_track','var')
        if TEST_probabilistic
            unisys_file_mat = strrep(strrep(unisys_file,'s.','s_'),'.txt','_prob.mat');
        else
            unisys_file_mat = strrep(strrep(unisys_file,'s.','s_'),'.txt','_hist.mat');
        end
        fprintf('loading tc tracks from %s\n',unisys_file_mat);
        load(unisys_file_mat)
    end
    
    % generate all the rain footprints and save hazard torrential rain set
    hazard = climada_tr_hazard_set(tc_track, hazard_set_file_tr, centroids);
    
    fprintf('TR: max(max(hazard.arr))=%f\n',full(max(max(hazard.arr)))); % a kind of easy check

else
    fprintf('loading hazard TR from %s\n',hazard_set_file_tr);
    load(hazard_set_file_tr)
end

% show biggest TR event
figure(main_fig);
subplot(3,1,2)
values                = full(hazard.arr(max_tc_pos,:)); % get one tc footprint
centroids.Longitude   = hazard.lon; % as the gridding routine needs centroids
centroids.Latitude    = hazard.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(values,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
plot(centroids.Longitude,centroids.Latitude,'.r','MarkerSize',1);
if isfield(centroids,'onLand')
    water_points=find(centroids.onLand==0);
    plot(centroids.Longitude(water_points),centroids.Latitude(water_points),'.b','MarkerSize',1);
end
box on
climada_plot_world_borders
axis equal
axis(centroids_rect);
colorbar
title('TC Rain Sum [mm]');
cmap = climada_colormap(hazard.peril_ID);
colormap(cmap)
if use_cbfreeze
    freezeColors %freeze this plot's colormap
    cbfreeze(colorbar)
end

if ~isempty(TEST_location)
    text(TEST_location.longitude,TEST_location.latitude,TEST_location.name)
    plot(TEST_location.longitude,TEST_location.latitude,'xk');
end


% 4) HAZARD TROPICAL CYCLONE SURGE
% =================================
% hazard on input: the tropical cyclone (TC) hazard event set
% hazard on output: the storm surge (TS) hazard event set
if ~exist(hazard_set_file_ts,'file') || force_recalc
    hazard = tc_surge_hazard_create(hazard_TC,hazard_set_file_ts);
else
    fprintf('loading hazard TS from %s\n',hazard_set_file_ts);
    load(hazard_set_file_ts)
end


% show biggest TS event
figure(main_fig);
subplot(3,1,3)
values                = full(hazard.arr(max_tc_pos,:)); % get one tc footprint
centroids.Longitude   = hazard.lon; % as the gridding routine needs centroids
centroids.Latitude    = hazard.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(values,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
plot(centroids.Longitude,centroids.Latitude,'.r','MarkerSize',1);
if isfield(centroids,'onLand')
    water_points=find(centroids.onLand==0);
    plot(centroids.Longitude(water_points),centroids.Latitude(water_points),'.b','MarkerSize',1);
end
box on
climada_plot_world_borders
axis equal
axis(centroids_rect);
colorbar
title('TC Surge Field [m]');
cmap = climada_colormap(hazard.peril_ID);
colormap(cmap)

if ~isempty(TEST_location)
    text(TEST_location.longitude,TEST_location.latitude,TEST_location.name)
    plot(TEST_location.longitude,TEST_location.latitude,'xk');
end

return
