function [hazard,EDS,centroids,entity]=tc_surge_Barisal(...
    adm_lvl,force_centroids_recalc,force_entity_recalc,force_hazard_recalc, check_plots)
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
%          call climada_ts_hazard_set in order to
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
% Gilles Stassen, gillesstassen@hotmail.com, 20150202, split entity and
% centroid generation into separate blocks; new input args: force_centroids_recalc, force_entity_recalc
% David N. Bresch, david.bresch@gmail.com, 20150819, centroids in their own dir
% Lea Mueller, muellele@gmail.com, 20151125, rename to climada_centroids_generate from climada_generate_centroids
% David N. Bresch, david.bresch@gmail.com, 20160529, calling climada_ts_hazard_set instead of tc_surge_hazard_create
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
if ~exist('force_entity_recalc','var'),             force_entity_recalc = 0;            end
if ~exist('force_centroids_recalc','var'),          force_centroids_recalc = 0;         end
if ~exist('force_hazard_recalc','var'),             force_hazard_recalc = 0;            end
if ~exist('check_plots','var'),                     check_plots = 1;                    end

% PARAMETERS
%
% set global variables (be careful, they should be reset, see bottom of code)
if ~isfield(climada_global,'climada_global_ori')
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
%
% The shape file with detailed border info
BGD_admin_regions_shapefile = [module_data_dir filesep 'entities' filesep ...
    'BGD_adm' filesep 'BGD_adm' num2str(adm_lvl) '.shp']; % adm_lvl is an input parameter
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
location.lon  = 90.3667;
location.lat  = 22.7000;
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

% 1) Centroids for study region
% -----------------------------
% Define the file with centroids (geo-locations of the points we later
% evaluate and store storm surge heights at), as well as the entity file
% directories.
% see climada_create_GDP_entity to create centroids file

centroids_file  =   [module_data_dir filesep 'centroids' filesep strcat('Barisal_',num2str(adm_lvl),'_centroids.mat')];
entity_file     =   [module_data_dir filesep 'entities' filesep strcat('Barisal_',num2str(adm_lvl),'_entity.mat')];
entity_file_xls =   [module_data_dir filesep 'entities' filesep strcat('Barisal_',num2str(adm_lvl),'_entity.xls')];
%hrnl_img_file   =   [module_data_dir filesep 'entities' filesep '89_21_91_24_Bangladesh_Barisal_F182010.v4c.avg_lights_x_pct.lzw.tiff'];


% 2) Tropical cyclone (TC) tracks
% -------------------------------
% Set UNISYS TC track data file (for info, see climada_tc_read_unisys_database)
unisys_file= [climada_global.data_dir filesep 'tc_tracks' filesep tc_track_file];
precip_data_file = [module_data_dir filesep 'precip_data' filesep 'precip.mon.total.v6.nc'];

% 3) Bathymetry parameters are set in climada_ts_hazard_set
srtm_data_dir = [module_data_dir filesep 'system' filesep 'srtm_55_08'];
DEM_save_file = [module_data_dir filesep 'system' filesep strcat('Barisal_',num2str(adm_lvl),'_DEM.mat')];

% 4) Surge hazard event set
% -------------------------
% Define the hazard event set file to store the Barisal TC and TS hazard
% event sets
hazard_set_file_tc=[module_data_dir filesep 'hazards' filesep strcat('Barisal_',num2str(adm_lvl),'_hazard_TC.mat')];
hazard_set_file_ts=[module_data_dir filesep 'hazards' filesep strcat('Barisal_',num2str(adm_lvl),'_hazard_TS.mat')];
hazard_set_file_tr=[module_data_dir filesep 'hazards' filesep strcat('Barisal_',num2str(adm_lvl),'_hazard_TR.mat')];
hazard_set_file_ma=[module_data_dir filesep 'hazards' filesep strcat('Barisal_',num2str(adm_lvl),'_hazard_MA.mat')];
hazard_set_file_mod_tr=[module_data_dir filesep 'hazards' filesep strcat('Barisal_',num2str(adm_lvl),'_hazard_mod_TR.mat')];
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
if isempty(force_centroids_recalc)
    force_centroids_recalc=0;
elseif force_centroids_recalc == 1 && force_hazard_recalc == 0
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
    if adm_lvl == 4
        centroids_90m_check = 1;
    else centroids_90m_check = 0; end
    if ~centroids_90m_check
        centroid_resolution_km = 1.8;
        centroids = climada_centroids_generate(centroids_rect,centroid_resolution_km);
        % Assign elevation to centroids
        [DEM, centroids] = climada_read_srtm_DEM('DL',centroids, DEM_save_file, 5, check_plots);
    else
        % generate centroids from DEM
        [DEM, centroids] = climada_read_srtm_DEM('DL',centroids_rect, DEM_save_file, 8, check_plots);
    end
    fprintf('saving centroids to %s \n',centroids_file)
    save(centroids_file, 'centroids')
end

if isempty(force_entity_recalc), force_entity_recalc=0; end

if (exist(entity_file,'file') || exist(entity_file_xls,'file')) && ~force_entity_recalc
    if exist(entity_file,'file') && exist(entity_file_xls,'file')
        xls_dir = dir(entity_file_xls); mat_dir = dir(entity_file);
        if xls_dir.datenum - mat_dir.datenum >0 % Read excel file only if it is newer
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
        entity_BGD = climada_nightlight_entity(country_name,'Barisal'); % 20160907, dnb
        save(entity_BGD_file,'entity_BGD');     % Save for next time
    end
    
    % Clip the centroids and entity to the bounding box (centroids_rect) of
    % the desired region. Increase the resolution of centroids (also
    % possible for assets, but not really necessary)
    entity = climada_entity_crop(entity_BGD, centroids_rect,1);
    
    % Encode each asset to nearest on-land centroid for damage calculations
    temp_centroids=centroids;
    temp_centroids.centroid_ID(centroids.onLand ~=1)= - temp_centroids.centroid_ID(centroids.onLand ~=1);
    entity.assets = climada_assets_encode(entity.assets,temp_centroids);
    clear temp_centroids % Free up memory
    
    save(entity_file, 'entity')
    climada_entity_save_xls(entity,entity_file_xls,1,0,0);
    % The last three input args define whether damage functions, measures
    % and discounts are overwritten, respectively. Assets are always
    % overwritten.
end


% Plot entity assets
if check_plots
    entity.assets.reference_year = 2014;
    figure('name','Asset distribution Barisal','color','w')
    climada_plot_entity_assets(entity);
    hold on
    % Plot centroids
    if isfield(centroids,'onLand')
        ndx = centroids.onLand ==1;
        plot(centroids.lon(ndx),centroids.lat(ndx),'.g','markersize',1);
        plot(centroids.lon(~ndx),centroids.lat(~ndx),'.b','markersize',1);
    end
    % Plot location
    % if ~isempty(location)
    %     text(location.longitude+0.05,location.latitude,10,location.name,'fontsize',14,'color',[1 0 0],'backgroundcolor',[1 1 1])
    %     plot3(location.longitude,location.latitude,10,'or','markersize',40, 'linewidth', 3);
    % end
    % Plot details
    if details_check
        climada_shp_explorer([module_data_dir filesep 'entities' filesep 'BGD_shapes' filesep 'buildings.shp']);
    end
    axis equal
    axis(centroids_rect)
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
            fprintf('No inland decay, consider module tc_hazard_advanced\n');
        end
        
        % Expand set of tracks by generating probabilistic tracks
        % See function header for more details on generating probabilistic
        % tracks, such as specifying ensemble size, max angle etc.
        tc_track = climada_tc_random_walk(tc_track);
        if exist('climada_tc_track_wind_decay_calculate','file')
            % Add the inland decay correction to all probabilistic nodes
            tc_track   = climada_tc_track_wind_decay(tc_track, p_rel,0);
        end
        
        % Plot the tracks
        if check_plots
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
    %save(hazard_set_file_tc,'hazard_tc');
    
else
    fprintf('loading TC wind hazard set from %s\n',hazard_set_file_tc);
    load(hazard_set_file_tc); % Load existing hazard
    if ~exist('hazard_tc','var')
        hazard_tc = hazard; clear hazard;
    end
end

% 3) Create TS hazard event set
% -----------------------------------
if ~exist(hazard_set_file_ts,'file') || force_hazard_recalc
    % policy research working paper 5280 "Vulnerability of Bangladesh to Cyclones in a Changing Climate"
    %surge_params = [0.1252 -1.7005]; % 20160529, no parames handed over to climada_ts_hazard_set
    hazard_ts = climada_ts_hazard_set(hazard_tc,hazard_set_file_ts);
else
    fprintf('loading TS surge hazard set from %s\n',hazard_set_file_ts);
    load(hazard_set_file_ts);
    if ~exist('hazard_ts','var')
        hazard_ts = hazard; clear hazard;
    end
end

% 4) Create TR hazard set
if ~exist(hazard_set_file_tr, 'file')  || force_hazard_recalc
    if ~exist('tc_track','var')
        tc_track = climada_tc_read_unisys_database(unisys_file);
        if probabilistic
            if exist('climada_tc_track_wind_decay_calculate','file')
                [~, p_rel]  = climada_tc_track_wind_decay_calculate(tc_track,0);
            else fprintf('No inland decay, consider module tc_hazard_advanced\n'); end
            
            tc_track = climada_tc_random_walk(tc_track);
            close
            if exist('climada_tc_track_wind_decay_calculate','file')
                tc_track   = climada_tc_track_wind_decay(tc_track, p_rel,0);
            end
        end
    end
    hazard_tr = climada_tr_hazard_set(tc_track, hazard_set_file_tr,centroids);
else
    fprintf('loading TR rain hazard set from %s\n',hazard_set_file_tr);
    load(hazard_set_file_tr);
    if ~exist('hazard_tr','var')
        hazard_tr = hazard; clear hazard;
    end
end

% 5) Validate TR hazard set using historical APHRODITE MA daily
%    precipitation data
if ~exist(hazard_set_file_ma,'file') || force_hazard_recalc
    hazard_ma = climada_ma_hazard_set(1972:2007,centroids,hazard_set_file_ma,check_plots);  
else
    fprintf('loading MA rain hazard set from %s \n',hazard_set_file_ma);
    load(hazard_set_file_ma);
    hazard_ma = hazard; clear hazard;
end
if ~exist(hazard_set_file_mod_tr,'file') || force_hazard_recalc
    hazard_mod_tr = climada_mod_tr_hazard_set(hazard_tr,hazard_ma,location,hazard_set_file_mod_tr, check_plots);   
else
    fprintf('loading modified TR hazard set file from %s \n',hazard_set_file_mod_tr);
    load(hazard_set_file_mod_tr);
    hazard_mod_tr = hazard; clear hazard;
end

if check_plots
    % Plot the hazard sets for largest event
    [~,max_tc_pos]=max(sum(hazard_tc.intensity,2).* hazard_tc.orig_event_flag'); % the maximum TC intensity
    max_event_fig_title = sprintf('Storm surge Barisal admin 3 for largest event # %i in the year %i',max_tc_pos,hazard_tc.yyyy(max_tc_pos));

    figure('Name',max_event_fig_title,'Color',[1 1 1])
    climada_hazard_plot_hr(hazard_tc,max_tc_pos);
    hold off
    
    figure('Name',max_event_fig_title,'Color',[1 1 1])
    climada_hazard_plot_hr(hazard_ts,max_tc_pos);
    hold off
    
    figure('Name',max_event_fig_title,'Color',[1 1 1])
    climada_hazard_plot_hr(hazard_tr,max_tc_pos);
    hold off
    drawnow
    
    % Plot the max intensity at each centroid
    figure('Name','Maximum hazard intenstity at each centroid','color','w')
    climada_hazard_plot_hr(hazard_tc,0);
    hold off
    
    figure('Name','Maximum hazard intenstity at each centroid','color','w')
    climada_hazard_plot_hr(hazard_ts,0);
    hold off
    
    figure('Name','Maximum hazard intenstity at each centroid','color','w')
    climada_hazard_plot_hr(hazard_tr,0);
    hold off
    drawnow
    
    % Plot intensity v. return period for storm surge: the integers specify the
    % return periods of interest
    climada_hazard_stats(hazard_ts,[1 4 10 40 100 400],1,'TS');
    hold off
end

% Generate damage set and adaptation cost curve
% ----------------------------------------------

% Ensure asset covers are set to asset values, deductibles set to zero for
% reasonable damage calculation (the equivalent of ignoring insurance
% policies altogether)
entity.assets.Cover = entity.assets.Value;
entity.assets.Deductible = entity.assets.Value .* 0;

entity = climada_assets_encode(entity,hazard_tc);

fprintf('calculating event damage sets for TC, TS and TR hazards \n')
EDS_tc = climada_EDS_calc(entity,hazard_tc);
EDS_ts = climada_EDS_calc(entity,hazard_ts);
%EDS_tr = climada_EDS_calc(entity,hazard_tr);

% tc_surge_plot_3d_Barisal(hazard_ts,max_tc_pos);

if check_plots
    climada_EDS_plot_3d(hazard_tc,EDS_tc);
    climada_EDS_plot_3d(hazard_ts,EDS_ts);
%    climada_plot_EDS_3d(hazard_tr,EDS_tr);
end
if check_plots
    fprintf('calculating impact of measures on event damage sets \n')
    TS_measures_impact = climada_measures_impact(entity,hazard_ts,'no');
    TS_measures_impact.title_str = 'Adaptation Cost Curve Barisal Province';
    climada_adaptation_cost_curve(TS_measures_impact,[],[],[],[],[],1); % The 1 triggers reverse benefit-cost ratio (as opposed to cost-benefit)
end

EDS.TS = EDS_ts;
EDS.TC = EDS_tc;
% EDS.TR = EDS_tr;

hazard.TS       = hazard_ts;
hazard.TC       = hazard_tc;
hazard.TR.ori   = hazard_tr;
hazard.TR.mod   = hazard_mod_tr;
hazard.MA       = hazard_ma;
% reset global variables
if isfield(climada_global, 'climada_global_ori')
    climada_global = climada_global.climada_global_ori;
end

return
