function [hazard,EDS,centroids,entity,profile_stats]=tc_surge_Barisal_HR(force_centroids_recalc,force_entity_recalc, force_hazard_recalc,check_plots)
% climada
% MODULE:
%   barisal_demo
% NAME:
%   tc_surge_Barisal_HR
% PURPOSE:
%   This code creates the tropical cyclone (TC) storm
%   surge (TS) hazard for an area around the city of Barisal, Bangladesh.
%       1) get centroids for the area of interest(see PARAMETERS)
%          if they do not exist, run GDP_entity in order to create them
%       2) create TC wind hazard event set
%          call tc_surge_hazard_create in order to
%       3) create bathymetry file for region
%     	4) create TC surge hazard event set
%   show the result
%
%   See tc_surge_plot_3d for 3D plots of surge fields in general
%
% CALLING SEQUENCE:
%   [hazard, EDS, centroids, entity] = tc_surge_Barisal
% EXAMPLE:
%   tc_surge_Barisal
%   [hazard, EDS, centroids, entity]=tc_surge_Barisal(1,0,1) - save output
%       to variables - forces the recalculation of hazard set, but not
%       entity/centroids.
% INPUTS:
% OPTIONAL INPUT PARAMETERS:
%   force_centroids_entity_recalc: Automatically set to 0. Set to 1 if you
%       wish to recalculate the centroids and entity, despite the relevant
%       files already existing - calculation will take longer.
%   force_hazard_recalc: Automatically set to 0. Set to 1 if you
%       wish to recalculate the TS and TC hazards, despite the relevant
%       files already existing - calculation will take longer.
% OUTPUTS:
%   hazard: Struct with fields
%             .TS:    Storm surge hazard set (with usual fields)
%             .TC:    Tropical cyclone hazard set
%             .TR:    .ori:   Torrential rain hazard set from rainfields 
%                             calculated usingR-CLIPER
%                     .mod:   Modified TR hazard set to fit with the
%                             historical precipitation information stored
%                             in the MA hazard set
%             .MA:    Monsoon Asia hazard set, containing only historical
%                     data from the APHRODITE data set - used for validation
%   EDS:    Struct with fields
%             .TS:    Storm surge event damage set
%             .TC:    Tropical cyclone event damage set
%   centroids:  High resolution centroids
%   entity:     High resolution entity
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20141121
% David N. Bresch, david.bresch@gmail.com, 20141215, cleanup, especially the map border stuff
% Gilles Stassen, 20150109, high resolution copy of tc_surge_Barisal,
%                           hardwired for admin 3
%-

hazard=[]; EDS = []; centroids = []; entity = []; profile_stats = [];% init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% Check input variables
if ~exist('force_centroids_recalc', 'var'), force_centroids_recalc = 0; end
if ~exist('force_entity_recalc',    'var'), force_entity_recalc = 0;    end
if ~exist('force_hazard_recalc',    'var'), force_hazard_recalc = 0;    end
if ~exist('check_plots',            'var'), check_plots = 1;            end
% PARAMETERS
%
% set global variables (be careful, they should be reset, see bottom of code)
if ~isfield(climada_global, 'climada_global_ori')
    climada_global.climada_global_ori = climada_global; % store for reset
end
climada_global.EDS_at_centroid = 1;
climada_global.waitbar = 0; % suppress waitbar
climada_global.tc.default_min_TimeStep = 1/4;
%
% the module's data folder:
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% in case one needs to access data of another module, use (eg for country_risk)
%country_risk_module_data_dir=[fileparts(fileparts(which('country_risk_calc'))) filesep 'data'];
%
% the shape file with higher resolution for Bangladesh
BGD_country_shapefile=[module_data_dir filesep 'entities' filesep 'BGD_adm' filesep 'BGD_adm0.shp'];

% the shape file with higher resolution borders/rivers/coast (choose one
% level up from study admin level as a buffer)
BGD_admin2_shapefile = [module_data_dir filesep 'entities' filesep ...
    'BGD_adm' filesep 'BGD_adm2.shp'];
%
% The shape file with detailed admin 4 info
BGD_admin3_shapefile = [module_data_dir filesep 'entities' filesep ...
    'BGD_adm' filesep 'BGD_adm3.shp'];
%
% If we want to add geographical details to the entity plot, set to 1
details_check = 0;
% name of the country
country_name    = 'Bangladesh';
%
% name of the tropical cyclone track file:
tc_track_file   = 'tracks.nio.txt';
%
% annotation (just to label Barisal on plots):
location.name   = '  Barisal'; % first two spaces for nicer labeling
location.longitude  = 90.3667;
location.latitude   = 22.7000;
%
% Define whether we run the simulation using only historical tracks, or
% generate the full probabilistic TC hazard (and subsequent TS hazard)
probabilistic = 1;

% make sure there is a .mat version of the border shapes file
[fP,fN] = fileparts(BGD_country_shapefile);
BGD_country_shapefile_mat=[fP filesep fN '.mat'];
if ~exist(BGD_country_shapefile_mat,'file')
    % first time, read the shape file and store as .mat
    climada_shaperead(BGD_country_shapefile);
end

% replace the map border shape file (be careful, they should be reset)
climada_global.map_border_file=BGD_country_shapefile_mat;

% % make sure there is a .mat version of the border shapes file
% [fP,fN] = fileparts(BGD_admin2_shapefile);
% BGD_admin2_shapefile_mat=[fP filesep fN '.mat'];
% if ~exist(BGD_admin2_shapefile_mat,'file')
%     % first time, read the shape file and store as .mat
%     climada_shaperead(BGD_admin2_shapefile);
% end

% % replace the map border shape file (be careful, they should be reset)
% orig_climada_global_map_border_file=climada_global.map_border_file; % store for reset
% climada_global.map_border_file=BGD_admin2_shapefile_mat;

profile clear; profile on

% 1) Centroids for study region
% -----------------------------
% Define the file with centroids (geo-locations of the points we later
% evaluate and store storm surge heights at), as well as the entity file
% directories.
% see climada_create_GDP_entity to create centroids file

centroids_file  =   [module_data_dir filesep 'system' filesep 'Barisal_HR_centroids.mat'];
entity_file     =   [module_data_dir filesep 'entities' filesep 'Barisal_HR_entity.mat'];
entity_file_xls =   [module_data_dir filesep 'entities' filesep 'Barisal_HR_entity.xls'];

% 2) Tropical cyclone (TC) tracks
% -------------------------------
% Set UNISYS TC track data file (for info, see climada_tc_read_unisys_database)
unisys_file= [climada_global.data_dir filesep 'tc_tracks' filesep tc_track_file];

% 3) Bathymetry parameters are set in tc_surge_hazard_create
srtm_data_dir = [module_data_dir filesep 'system' filesep 'srtm_55_08'];
DEM_save_file = [module_data_dir filesep 'system' filesep 'Barisal_HR_DEM.mat'];

% 4) Surge hazard event set
% -------------------------
% Define the hazard event set file to store the Barisal TC and TS hazard
% event sets
hazard_set_file_tc=[module_data_dir filesep 'hazards' filesep 'Barisal_HR_hazard_TC.mat'];
hazard_set_file_ts=[module_data_dir filesep 'hazards' filesep 'Barisal_HR_hazard_TS.mat'];
hazard_set_file_tr=[module_data_dir filesep 'hazards' filesep 'Barisal_HR_hazard_TR.mat'];
hazard_set_file_ma=[module_data_dir filesep 'hazards' filesep 'Barisal_HR_hazard_MA.mat'];
hazard_set_file_mod_tr=[module_data_dir filesep 'hazards' filesep 'Barisal_HR_hazard_mod_TR.mat'];
% CALCULATIONS
% ==============
% 1) Read the centroids
% ---------------------

% If the file specified above exists, use the high resolution borders
% stored there, else revert to default climada low resolution borders.
if exist(BGD_admin3_shapefile,'file')
    BGD_admin_regions = shaperead(BGD_admin3_shapefile);
    
    for i = 1 : numel(BGD_admin_regions)
        if strfind(BGD_admin_regions(i).NAME_3, 'Barisal')
            ID = BGD_admin_regions(i).ID_3;
            break;
        end
        if i == numel(BGD_admin_regions), fprintf('ERROR: Barisal not found'); return; end
    end
    
    % Prep the region we need from the bounding box
    centroids_rect =[BGD_admin_regions(i).BoundingBox(:,1)' BGD_admin_regions(i).BoundingBox(:,2)'];
else
    fprintf('WARNING: Bangladesh shapefiles not found, using standard resolution world map instead\n');
    fprintf('Download shapefiles from: www.diva-gis.org/gdata \n');
    if exist('shapes', 'var')
        centroids_rect =[shapes.BoundingBox(:,1)' shapes.BoundingBox(:,2)'];
    end
end

% Load existing centroids and entity files if it exists, and unless user
% specifies to recaulculate the centroids & entity
if isempty(force_centroids_recalc)
    force_centroids_recalc=0;
elseif force_centroids_recalc && ~force_hazard_recalc;
    cprintf([0.25 0.25 1],'NOTE: recalculation of centroids mandates regeneration of hazard sets \n');
    force_hazard_recalc = 1;
end

if exist(centroids_file,'file') && exist(DEM_save_file,'file') &&  ~force_centroids_recalc
    fprintf('loading centroids from %s \n',centroids_file)
    load(centroids_file)    % load centroids
    fprintf('loading digital elevation model from %s \n', DEM_save_file)
    load(DEM_save_file)     % load DEM
else
    % Get centroids for the whole of Bangladesh
    % centroids_BGD = climada_create_GDP_entity(country_name,[],0,1);
    
    % Get high resolution centroids for bounding box
    centroids_90m_check = 1;
    if ~centroids_90m_check
        % Generate centroids independently of DEM (only sensible when
        % choosing a lower resolution than that of DEM)
        centroids = climada_generate_HR_centroids(centroids_rect,1.8);
        % assign elevation to centroids
        [DEM, centroids] = climada_read_srtm_DEM(srtm_data_dir,centroids, DEM_save_file, 5, 0);
    else        
        [DEM, centroids] = climada_read_srtm_DEM(srtm_data_dir,centroids_rect, DEM_save_file, 5, 0);
    end
    fprintf('saving centroids to %s \n',centroids_file)
    save(centroids_file, 'centroids')
end

if isempty(force_entity_recalc), force_entity_recalc=0; end

if (exist(entity_file,'file') || exist(entity_file_xls,'file')) && ~force_entity_recalc
    if exist(entity_file,'file') && exist(entity_file_xls,'file')
        xls_dir = dir(entity_file_xls); mat_dir = dir(entity_file);
        if datenum(xls_dir.date) > datenum(mat_dir.date) % Read excel file only if it is newer
            fprintf('loading entity from %s \n',entity_file_xls)
            entity = climada_entity_read(entity_file_xls,'SKIP');  % Read existing entity Excel file
            save(entity_file,'entity') % Save any updates made to excel file in mat file
        else
            fprintf('loading entity from %s \n',entity_file)
            load(entity_file)
        end
    else
        fprintf('loading entity from %s \n',entity_file)
        load(entity_file)
    end % load entity
else
    % Get enitity for the whole of Bangladesh (such that the asset values
    % are scales to GDP
    entity_BGD_file = [module_data_dir filesep 'entities' filesep 'BGD_HR_entity.mat'];
    if exist(entity_BGD_file,'file')
        load(entity_BGD_file)       % Load existing file
    else
        % Entity from high resolution night lights
        entity_BGD = climada_nightlight_entity(country_name,'Barisal',1);
        save(entity_BGD_file,'entity_BGD');     % Save for next time
    end
    
    % Clip the centroids and entity to the bounding box (centroids_rect) of
    % the desired region. Increase the resolution of centroids (also
    % possible for assets, but not really necessary)
    % Clip and increase resolution of centroids. Waterborne assets moved onto land
    [~, entity] = climada_clip_centroids_entity([], entity_BGD, centroids_rect, [], 1);
    
    % Encode each asset to nearest on-land centroid for damage calculations
    temp_centroids=centroids;
    temp_centroids.centroid_ID(centroids.onLand ~=1)= - temp_centroids.centroid_ID(centroids.onLand ~=1);
    entity.assets = climada_assets_encode(entity.assets,temp_centroids);
    clear temp_centroids % Free up memory
    
    save(entity_file, 'entity')
    climada_entity_save_xls(entity,entity_file_xls,0,0,0);
    % The last three input args define whether damage functions, measures
    % and discounts are overwritten, respectively. Assets are always
    % overwritten.
end

% Plot entity assets
if check_plots
    entity.assets.reference_year = 2014;
    figure('name','Asset distribution Barisal','color','w')
    colormap(flipud(copper))
    climada_plot_entity_assets(entity,centroids);
    hold on
    % Plot location
%     if ~isempty(location)
%         text(location.longitude+0.05,location.latitude,10,location.name,'fontsize',14,'color',[1 0 0],'backgroundcolor',[1 1 1])
%         plot3(location.longitude,location.latitude,10,'or','markersize',40, 'linewidth', 3);
%     end
    % Plot details
    if details_check
        climada_shp_explorer([module_data_dir filesep 'entities' filesep 'BGD_shapes' filesep 'buildings.shp']);
    end
    hold off
end

% 2) Create TC hazard event set
% -----------------------------------
% Do complete calculation if no TC hazard set file exists, or if demanded
% by user, else load existing file.
if ~exist(hazard_set_file_tc,'file') || force_hazard_recalc
    % Load historical tracks
    tc_track = climada_tc_read_unisys_database(unisys_file);
    if probabilistic % Do complete probabilistic calculation
        if exist('climada_tc_track_wind_decay_calculate','file')
            % Wind speed decay at track nodes after landfall
            [~, p_rel]  = climada_tc_track_wind_decay_calculate(tc_track,0);
        else
            fprintf('NO inland decay, consider module tc_hazard_advanced\n');
        end
        
        % Expand set of tracks by generating probabilistic tracks
        % See function header for more details on generating probabilistic
        % tracks, such as specifying ensemble size, max angle etc.
        tc_track = climada_tc_random_walk_position_windspeed(tc_track,[],[],[],[],0);
        close
        if exist('climada_tc_track_wind_decay_calculate','file')
            % Add the inland decay correction to all probabilistic nodes
            tc_track   = climada_tc_track_wind_decay(tc_track, p_rel,0);
        end
        if check_plots
            % Plot the tracks
            figure('Name','TC tracks','Color',[1 1 1]);
            hold on
            for event_i=1:length(tc_track) % plot all tracks
                plot(tc_track(event_i).lon,tc_track(event_i).lat,'-b');
            end % event_i
            % Overlay historic (to make them visible, too)
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
        end
    end
    % Generate all the wind footprints: create TC hazard set
    hazard_tc = climada_tc_hazard_set(tc_track, hazard_set_file_tc, centroids);
    hazard_tc.units = 'm/s'; % Set the units for the plot
    save(hazard_set_file_tc,'hazard_tc');
    
else
    fprintf('loading TC wind hazard set from %s \n',hazard_set_file_tc);
    load(hazard_set_file_tc); % Load existing hazard        
    %hazard_tc = hazard; clear hazard;
end

% 3) Create TS hazard event set
% -----------------------------------
if ~exist(hazard_set_file_ts,'file') || force_hazard_recalc
    hazard_ts = tc_surge_hazard_create(hazard_tc,hazard_set_file_ts,DEM);
else
    fprintf('loading TS hazard set file from %s \n',hazard_set_file_ts);
    load(hazard_set_file_ts);
    hazard_ts = hazard; clear hazard;
end

% 4) Create TR hazard set
if ~exist(hazard_set_file_tr, 'file')  || force_hazard_recalc
    if ~exist('tc_track','var')
        tc_track = climada_tc_read_unisys_database(unisys_file);
        if probabilistic
            if exist('climada_tc_track_wind_decay_calculate','file')
                [~, p_rel]  = climada_tc_track_wind_decay_calculate(tc_track,0);
            else fprintf('No inland decay, consider module tc_hazard_advanced\n'); end
            
            tc_track = climada_tc_random_walk_position_windspeed(tc_track,[],[],[],[],0,0);
            close
            if exist('climada_tc_track_wind_decay_calculate','file')
                tc_track   = climada_tc_track_wind_decay(tc_track, p_rel,0);
            end
        end
    end
    hazard_tr = climada_tr_hazard_set(tc_track, hazard_set_file_tr,centroids);
else
    fprintf('loading TR hazard set file from %s \n',hazard_set_file_tr);
    load(hazard_set_file_tr);
    hazard_tr = hazard; clear hazard;
end

% 5) Validate TR hazard set using historical APHRODITE MA daily
%    precipitation data
if ~exist(hazard_set_file_ma,'file') || force_hazard_recalc
    hazard_ma = climada_ma_hazard_set(1972:2007,centroids,hazard_set_file_ma,check_plots);  
else
    fprintf('loading MA hazard set file from %s \n',hazard_set_file_ma);
    load(hazard_set_file_ma);
    hazard_ma = hazard; clear hazard;
end
if ~exist(hazard_set_file_mod_tr,'file') || force_hazard_recalc
    hazard_mod_tr = climada_mod_tr_hazard_set(hazard_tr,hazard_ma,location,hazard_set_file_mod_tr, check_plots);   
else
    fprintf('loading RF hazard set file from %s \n',hazard_set_file_mod_tr);
    load(hazard_set_file_mod_tr);
    hazard_mod_tr = hazard; clear hazard;
end

if check_plots
    % The following maps define the colour spectra, downloaded from https://www.ncl.ucar.edu/Document/Graphics/ColorTables/MeteoSwiss.shtml
    % ts colours
    map_ts = [
        254 254 254
        223 255 249
        154 217 202
        103 194 163
        64 173 117
        50 166 150
        90 160 205
        66 146 199
        76 141 196
        7  47 107
        7  30  70
        76   0 115]./255;
    % tc colours
    map_tc = [
        255 255 255
        239 244 209
        232 244 158
        170 206  99
        226 237  22
        255 237   0
        255 237 130
        244 209 127
        237 165  73
        229 140  61
        219 124  61
        239   7  61
        232  86 163
        155 112 168
        99 112 247
        127 150 255
        142 178 255
        181 201 255]./255;
    % map_tc = [
    %     255 255 255
    %     255 245 204
    %     255 230 112
    %     255 204  51
    %     255 175  51
    %     255 153  51
    %     255 111  51
    %     255  85   0
    %     230  40  30
    %     200  30  20]./255;
    % tr colours
    map_tr = [linspace(1,0,12)' linspace(1,0,12)' linspace(1,0.75,12)'];
    
    % Plot the hazard sets for largest event
    [~,max_tc_pos]=max(sum(hazard_tc.intensity,2).* hazard_tc.orig_event_flag'); % the maximum TC intensity
    max_event_fig_title = sprintf('Hazard Barisal admin 3 for largest TC event # %i in the year %i',max_tc_pos,hazard_tc.yyyy(max_tc_pos));
    %main_fig=figure('Name',max_event_fig_title,'Position',[89 223 1521 413],'Color',[1 1 1]); % [89 223 1014 413]
    figure('Name',max_event_fig_title, 'color', 'w')
    %subplot(1,3,1)%'title', 'Max Wind Speed [m/s]')
    %h(1) = subaxis(1,3,1,'SV',0.1);
    colormap(map_tc) % set the colormap to match the max event plot
    climada_hazard_plot(hazard_tc,max_tc_pos,location);
    %cbfreeze
    %freezeColors(h(1))
    %figure(main_fig)
    figure('Name',max_event_fig_title, 'color', 'w')
    %subplot(1,3,2)%,'title', 'Max Storm Surge Height [m]')
    %h(2) = subaxis(1,3,2,'SV',0.1);
    colormap(map_ts) % Set the colormap to match the max event plot
    climada_hazard_plot(hazard_ts,max_tc_pos,location);
    %cbfreeze
    %freezeColors(h(2))
    %figure(main_fig)
    figure('Name',max_event_fig_title, 'color', 'w')
    %subplot(1,3,3)
    %h(3) = subaxis(1,3,3,'SV',0.1);
    colormap(map_tr) % Set the colormap to match the max event plot
    climada_hazard_plot(hazard_tr,max_tc_pos,location);
    hold off
    
    
    % Plot the max intensity at each centroid
    %max_intensity_fig = figure('Name','Maximum hazard intenstity at each centroid','Position',[89 223 1521 413],'Color',[1 1 1]); % [89 223 1014 413]
    %figure(max_intensity_fig)
    figure('Name','Maximum hazard intenstity at each centroid','color','w')
    %subplot(1,3,1)%'title', 'Max Wind Speed [m/s]')
    colormap(map_tc) % set the colormap to match the max event plot
    climada_hazard_plot(hazard_tc,0,location);
    %cbfreeze
    %freezeColors
    %figure(max_intensity_fig)
    figure('Name','Maximum hazard intenstity at each centroid','color','w')
    %subplot(1,3,2)%,'title', 'Max Storm Surge Height [m]')
    colormap(map_ts) % Set the colormap to match the max event plot
    climada_hazard_plot(hazard_ts,0,location);
    %cbfreeze
    %freezeColors
    %figure(max_intensity_fig)
    figure('Name','Maximum hazard intenstity at each centroid','color','w')
    %subplot(1,3,3)
    colormap(map_tr) % Set the colormap to match the max event plot
    climada_hazard_plot(hazard_tr,0,location);
    hold off
    
    
    % Plot intensity v. return period for storm surge: the integers specify the
    % return periods of interest
    climada_hazard_stats(hazard_ts,[1 4 10 40 100 400],1,'TS');
    hold off
    climada_hazard_stats(hazard_tc,[1 4 10 40 100 400],1,'TC');
    hold off
    climada_hazard_stats(hazard_tr,[1 4 10 40 100 400],1,'TR');
    hold off
end

% Generate damage set and adaptation cost curve
% ----------------------------------------------

% Ensure asset covers are set to asset values, deductibles set to zero for
% reasonable damage calculation (the equivalent of ignoring insurance
% policies altogether)
entity.assets.Cover = entity.assets.Value;
entity.assets.Deductible = entity.assets.Value .* 0;

EDS_tc = climada_EDS_calc(entity,hazard_tc);
EDS_ts = climada_EDS_calc(entity,hazard_ts);
EDS_tr = climada_EDS_calc(entity,hazard_tr);

% tc_surge_plot_3d_Barisal(hazard_ts,max_tc_pos);

if check_plots
    % 3d plot of event damage for largest event
    climada_plot_EDS_3d(hazard_tc,EDS_tc);
    climada_plot_EDS_3d(hazard_ts,EDS_ts);
    climada_plot_EDS_3d(hazard_tr,EDS_tr);
    
    % Cost benefit analysis of adaptation measures, and c-b plot
    TS_measures_impact = climada_measures_impact(entity,hazard_ts,'no');
    TS_measures_impact.title_str = 'Adaptation Cost Curve Barisal Province';
    climada_adaptation_cost_curve(TS_measures_impact,[],[],[],[],[],1); % The 1 triggers reverse benefit-cost ratio (as opposed to cost-benefit)
end

profile_stats = profile('INFO'); profile off;

EDS.TS = EDS_ts;
EDS.TC = EDS_tc;
EDS.TR = EDS_tr;

hazard.TS       = hazard_ts;
hazard.TC       = hazard_tc;
hazard.TR.ori   = hazard_tr;
hazard.TR.mod   = hazard_mod_tr;
hazard.MA       = hazard_ma;

% reset global variables
if isfield(climada_global,'climada_global_ori')
    climada_global = climada_global.climada_global_ori;
end


return
