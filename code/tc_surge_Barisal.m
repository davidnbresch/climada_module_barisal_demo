function [hazard,EDS,centroids,entity]=tc_surge_Barisal(adm_lvl,force_centroids_entity_recalc,force_hazard_recalc)
% climada
% MODULE:
%   barisal_demo
% NAME:
%   tc_surge_Barisal
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
%   adm_lvl: Specify the admin level of interest
%       Default:set to level 2
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
%   EDS:    Struct with fields
%             .TS:    Storm surge event damage set
%             .TC:    Tropical cyclone event damage set
%   centroids:  High resolution centroids
%   entity:     High resolution entity
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20141121
% David N. Bresch, david.bresch@gmail.com, 20141215, cleanup, especially the map border stuff
%-

hazard=[]; EDS = []; centroids = []; entity = []; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% Check input variables
if ~exist('adm_lvl','var'), adm_lvl = 2; end
if adm_lvl > 4
    fprintf('ERROR: Please choose an admin level from 0 to 4, corresponding to:\n');
    fprintf('\t 0: Country \n \t 1: Division \n \t 2: Province \n \t 3: District \n \t 4: Sub-district \n');
    return;
end
if ~exist('force_centroids_entity_recalc','var'), force_centroids_entity_recalc = 0; end
if ~exist('force_hazard_recalc','var'), force_hazard_recalc = 0; end

% PARAMETERS
%
% set global variables (be careful, they should be reset, see bottom of code)
climada_global_EDS_at_centroid=climada_global.EDS_at_centroid; % store for reset
climada_global.EDS_at_centroid = 1;
climada_global_waitbar = climada_global.waitbar; % store for reset
climada_global.waitbar = 0; % suppress waitbar
%
% the module's data folder:
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
%
% in case one needs to access data of another module, use (eg for country_risk)
%country_risk_module_data_dir=[fileparts(fileparts(which('country_risk_calc'))) filesep 'data'];
%
% the shape file with higher resolution for Bangladesh
BGD_country_shapefile=[module_data_dir filesep 'entities' filesep 'BGD_adm' filesep 'BGD_adm0.shp'];
%
% The shape file with detailed border info
BGD_admin_regions_shapefile = [module_data_dir filesep 'entities' filesep ...
    'BGD_adm' filesep 'BGD_adm' num2str(adm_lvl) '.shp']; % adm_lvl is an input parameter
%
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
orig_climada_global_map_border_file=climada_global.map_border_file; % store for reset
climada_global.map_border_file=BGD_country_shapefile_mat;

% 1) Centroids for study region
% -----------------------------
% Define the file with centroids (geo-locations of the points we later
% evaluate and store storm surge heights at), as well as the entity file
% directories.
% see climada_create_GDP_entity to create centroids file

centroids_file  =   [module_data_dir filesep 'system' filesep strcat('Barisal_',num2str(adm_lvl),'_centroids.mat')];
entity_file     =   [module_data_dir filesep 'entities' filesep strcat('Barisal_',num2str(adm_lvl),'_entity.mat')];
entity_file_xls =   [module_data_dir filesep 'entities' filesep strcat('Barisal_',num2str(adm_lvl),'_entity.xls')];
%hrnl_img_file   =   [module_data_dir filesep 'entities' filesep '89_21_91_24_Bangladesh_Barisal_F182010.v4c.avg_lights_x_pct.lzw.tiff'];


% 2) Tropical cyclone (TC) tracks
% -------------------------------
% Set UNISYS TC track data file (for info, see climada_tc_read_unisys_database)
unisys_file= [climada_global.data_dir filesep 'tc_tracks' filesep tc_track_file];

% 3) Bathymetry parameters are set in tc_surge_hazard_create

% 4) Surge hazard event set
% -------------------------
% Define the hazard event set file to store the Barisal TC and TS hazard
% event sets
hazard_set_file_tc=[module_data_dir filesep 'hazards' filesep strcat('Barisal_',num2str(adm_lvl),'_hazard_TC.mat')];
hazard_set_file_ts=[module_data_dir filesep 'hazards' filesep strcat('Barisal_',num2str(adm_lvl),'_hazard_TS.mat')];

% CALCULATIONS
% ==============
% 1) Read the centroids
% ---------------------
% Specify the administrative level of interest (from function input);
%       0 - Country
%       1 - Division
%       2 - Province
%       3 - District
%       4 - Sub-district
% The data was originally downloaded from http://www.diva-gis.org/gdata

% If the file specified above exists, use the high resolution borders
% stored there, else revert to default climada low resolution borders.
if exist(BGD_admin_regions_shapefile,'file')
    BGD_admin_regions = shaperead(BGD_admin_regions_shapefile);
    
    if adm_lvl > 0
        for i = 1 : numel(BGD_admin_regions)
            if strfind(eval(strcat('BGD_admin_regions(i).NAME_',num2str(adm_lvl))), 'Barisal')
                ID = eval(strcat('BGD_admin_regions(i).ID_',num2str(adm_lvl)));
                BGD_bounding_box = BGD_admin_regions(i).BoundingBox;
                break;
            end
            if i == numel(BGD_admin_regions), fprintf('ERROR: Barisal not found'); return; end
        end
    else
        BGD_bounding_box = BGD_admin_regions.BoundingBox;
    end
    
    % Prep the region we need from the bounding box
    centroids_rect =[BGD_bounding_box(:,1)' BGD_bounding_box(:,2)'];
else
    fprintf('WARNING: Bangladesh shapefiles not found, using standard resolution world map instead\n');
    fprintf('Download shapefiles from: www.diva-gis.org/gdata \n');
end

% Load existing centroids and entity files if it exists, and unless user
% specifies to recaulculate the centroids & entity
if isempty(force_centroids_entity_recalc),force_centroids_entity_recalc=0;end
if exist(centroids_file,'file') && ...
        (exist(entity_file,'file') || exist(entity_file_xls,'file')) && ...
        ~force_centroids_entity_recalc
    load(centroids_file)    % load centroids
    
    if exist(entity_file,'file') && exist(entity_file_xls,'file')
        xls_dir = dir(entity_file_xls); mat_dir = dir(entity_file);
        if datenum(xls_dir.date) > datenum(mat_dir.date) % Read excel file only if it is newer
            entity = climada_entity_read(entity_file_xls,'SKIP');  % Read existing entity Excel file
            % save(entity_file,'entity') % Save any updates made to excel file in mat file
        else
            load(entity_file)
        end
    else
        load(entity_file)
    end % load entity
else
    % Get enitity for the whole of Bangladesh (such that the asset values
    % are scales to GDP
    entity_BGD_file = [module_data_dir filesep 'entities\BGD_HR_entity.mat'];
    if exist(entity_BGD_file,'file')
        load(entity_BGD_file)       % Load existing file
    else
        % Entity from high resolution night lights
        entity_BGD = climada_high_res_entity(country_name,'Barisal',1);
        save(entity_BGD_file,'entity_BGD');     % Save for next time
    end
    
    % Get centroids for the whole of Bangladesh
    centroids_BGD = climada_create_GDP_entity(country_name,[],0,1);
    
    % Clip the centroids and entity to the bounding box (centroids_rect) of
    % the desired region. Increase the resolution of centroids (also
    % possible for assets, but not really necessary)
    fprintf('Clipping and upscaling resolution of centroids and entity assets... ');
    c_scale = [2 4 8 8 16]; % set scaling factor for centroids dependent on chosen admin level
    % Clip aand increase resolution of centroids. Waterborne assets moved onto land
    [centroids, entity] = climada_clip_centroids_entity(...
        centroids_BGD, entity_BGD, centroids_rect, c_scale(adm_lvl +1), 1);
    fprintf('Done \n');
    
    % Encode each asset to nearest on-land centroid for damage calculations
    entity.assets = climada_assets_encode_centroids(entity.assets,centroids,1,'centroids.onLand ==1');
    
    % Save centroids as .mat, entity as .mat and .xls
    save(centroids_file, 'centroids')
    save(entity_file, 'entity')
    climada_entity_save_xls(entity,entity_file_xls,0,0,0);
    % The last three input args define whether damage functions, measures
    % and discounts are overwritten, respectively. Assets are always
    % overwritten.
end

% Plot entity assets
entity.assets.reference_year = 2014;
figure('name','Asset distribution Barisal','color','w')
climada_plot_entity_assets(entity,centroids);
hold on
% Plot location
if ~isempty(location)
    text(location.longitude+0.05,location.latitude,10,location.name,'fontsize',14,'color',[1 0 0],'backgroundcolor',[1 1 1])
    plot3(location.longitude,location.latitude,10,'or','markersize',40, 'linewidth', 3);
    hold off
end

% 2) Create TC hazard event set
% -----------------------------------
if isempty(force_hazard_recalc),force_hazard_recalc = 0;end
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
        tc_track = climada_tc_random_walk(tc_track);
        close
        if exist('climada_tc_track_wind_decay_calculate','file')
            % Add the inland decay correction to all probabilistic nodes
            tc_track   = climada_tc_track_wind_decay(tc_track, p_rel,0);
        end
        
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
    % Generate all the wind footprints: create TC hazard set
    hazard_tc = climada_tc_hazard_set(tc_track, hazard_set_file_tc, centroids);
    hazard_tc.units = 'm/s'; % Set the units for the plot
    save(hazard_set_file_tc,'hazard_tc');
    
else
    fprintf('loading TC wind hazard set from %s\n',hazard_set_file_tc);
    load(hazard_set_file_tc); % Load existing hazard
end

fprintf('TC: max(max(hazard.intensity))=%f\n',full(max(max(hazard_tc.intensity)))); % a kind of easy check

% Show biggest TC event
[~,max_tc_pos]=max(sum(hazard_tc.intensity,2)); % the maximum TC intensity
main_fig=figure('Name',strcat('Storm surge Barisal admin ',num2str(adm_lvl)),'Position',[89 223 1014 413],'Color',[1 1 1]);
subplot(1,2,1)
values=full(hazard_tc.intensity(max_tc_pos,:)); % get one TC footprint
centroids.Longitude=hazard_tc.lon; % as the gridding routine needs centroids
centroids.Latitude=hazard_tc.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(values,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
% Uncomment the following 4 lines should you wish to plot centroid points
% plot(centroids.Longitude,centroids.Latitude,'.r','MarkerSize',1);
% if isfield(centroids,'onLand')
%     plot(centroids.Longitude(centroids.onLand == 0),centroids.Latitude(centroids.onLand == 0),'.b','MarkerSize',1);
% end
box on
climada_plot_world_borders;

axis equal
axis(centroids_rect);
% The following map defines the colour spectrum, downloaded from https://www.ncl.ucar.edu/Document/Graphics/ColorTables/MeteoSwiss.shtml
map_tc = [
    255 255 255
    255 245 204
    255 230 112
    255 204  51
    255 175  51
    255 153  51
    255 111  51
    255  85   0
    230  40  30
    200  30  20
    ]./255;
colormap(map_tc)
% We set a different colormap for the surge hazard plot later on, so freeze
% the colors in the existing sub plot.
freezeColors
cbfreeze(colorbar)
title(sprintf('Windfield [m/s] (event %i)',max_tc_pos));
fprintf('Max event %i\n',max_tc_pos);

% Plot location
if ~isempty(location)
    text(location.longitude,location.latitude,location.name)
    plot(location.longitude,location.latitude,'xk');
end

% 3) Create TS hazard event set
% -----------------------------------
if ~exist(hazard_set_file_ts,'file')
    hazard_ts = tc_surge_hazard_create(hazard_tc,hazard_set_file_ts,0,0);
else
    load(hazard_set_file_ts);
    hazard_ts = hazard; clear hazard;
end
% show biggest TS event
figure(main_fig);
subplot(1,2,2)
values = full(hazard_ts.intensity(max_tc_pos,:)); % get one tc footprint
centroids.Longitude=hazard_ts.lon; % as the gridding routine needs centroids
centroids.Latitude=hazard_ts.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(values,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
% Uncomment the following 4 lines should you wish to plot centroid points
% plot(centroids.Longitude,centroids.Latitude,'.r','MarkerSize',1);
% if isfield(centroids,'onLand')
%     plot(centroids.Longitude(centroids.onLand==0),centroids.Latitude(centroids.onLand==0),'.b','MarkerSize',1);
% end
box on
climada_plot_world_borders

axis equal
axis(centroids_rect);
% The following map defines the colour spectrum, downloaded from https://www.ncl.ucar.edu/Document/Graphics/ColorTables/MeteoSwiss.shtml
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
    76   0 115
    ]./255;
colormap(map_ts);
colorbar
title('Surgefield [m]');

% Plot location
if ~isempty(location)
    text(location.longitude,location.latitude,location.name)
    plot(location.longitude,location.latitude,'xk');
end

hold off

% Plot intensity v. return period: the integers specify the return periods
% of interest
climada_hazard_stats(hazard_ts,[1 4 10 40 100 400],1,'TS');

% Plot the max intensity at each centroid
max_intensity_fig = figure('Name','Maximum hazard intenstity at each centroid','Position',[89 223 1014 413],'Color',[1 1 1]);
figure(max_intensity_fig)
subplot(1,2,1)%'title', 'Max Wind Speed [m/s]')
colormap(map_tc) % set the colormap to match the max event plot
climada_hazard_plot(hazard_tc,0,location);
cbfreeze
freezeColors
figure(max_intensity_fig)
subplot(1,2,2)%,'title', 'Max Sturm Surge Height [m]')
colormap(map_ts) % Set the colormap to match the max event plot
climada_hazard_plot(hazard_ts,0,location);
hold off

% Generate damage set and adaptation cost curve
% ----------------------------------------------

% Ensure asset covers are set to asset values, deductibles set to zero for
% reasonable damage calculation (the equivalent of ignoring insurance
% policies altogether)
entity.assets.Cover = entity.assets.Value;
entity.assets.Deductible = entity.assets.Value .* 0;

% EDS_tc = climada_EDS_calc(entity,hazard_tc);
EDS_ts = climada_EDS_calc(entity,hazard_ts);

% tc_surge_plot_3d_Barisal(hazard_ts,max_tc_pos);

% climada_plot_EDS_3d(hazard_tc,EDS_tc);
climada_plot_EDS_3d(hazard_ts,EDS_ts);

TS_measures_impact = climada_measures_impact(entity,hazard_ts,'no');
TS_measures_impact.title_str = 'Adaptation Cost Curve Barisal Province';
climada_adaptation_cost_curve(TS_measures_impact,[],[],[],[],[],1); % The 1 triggers reverse benefit-cost ratio (as opposed to cost-benefit)

EDS.TS = EDS_ts;
% EDS.TC = EDS_tc;

hazard.TS = hazard_ts;
% hazard.TC = hazard_tc;

% reset global variables
climada_global.EDS_at_centroid=climada_global_EDS_at_centroid;
climada_global.waitbar=climada_global_waitbar;
climada_global.map_border_file=orig_climada_global_map_border_file;

return
