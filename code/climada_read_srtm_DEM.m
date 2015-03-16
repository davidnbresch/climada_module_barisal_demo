function [DEM, centroids] = climada_read_srtm_DEM(srtm_dir, centroidsORcountry, DEM_save_file, smooth, check_plot)
% climada
% MODULE:
%   barisal_demo
% NAME:
%   climada_read_srtm_DEM
% PURPOSE:
%   Read the digital elevation model data from the files in an existing
%   srtm directory. Data can be downloaded from http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp
% CALLING SEQUENCE:
%   DEM = climada_read_srtm_DEM(srtm_dir, centroids, DEM_save_file, smooth, check_plot)
% EXAMPLE:
%   DEM = climada_read_srtm_DEM(srtm_dir,[],[],[],1)
%   DEM = climada_read_srtm_DEM
%   [DEM, centroids] = climada_read_srtm_DEM(srtm_dir,[min_lon max_lon min_lat max_lat], DEM_save_file, 4,1)
% INPUTS:
%   srtm_dir:   The directory of an srtm data tile folder, containing at
%               least a .hdr and a .tif file. Can also be set to 'DL' which
%               will initiate automatic download from SRTM website
%               according to the centroid_rect given as input.
% OPTIONAL INPUT PARAMETERS:
%   centroids:  If centroids are provided as an input, the DEM will contain
%               elevation data sampled at the location of centroids, and
%               hence will have an extra field .centroid_ID. The extra
%               field is required for tc_surge_hazard_create if you wish to
%               provide your own topography data (e.g. srtm).
%               If this input is left empty , a
%               centroids struct will be generated at the same resolution
%               as the DEM.
%               If set to a 4-element vector (centroids_rect), these 4
%               points will define the area of interest, which is
%               subsequently cropped out of the DEM.
%               NOTE:   It is only sensible to provide a centroids struct
%               as input if its resolution is significantly lower than that
%               of the DEM, otherwise, it is much faster to generate
%               centroids directly from the DEM.
%   smooth:     Can either be set to an integer N (smooth by default filter
%               specified by a matrix size NxN with values 1/N^2) or a
%               smoothing filter. Default = [] (no smoothing).
%   check plot: Specify whether to plot a relief of the DEM, default = 0
% OUTPUTS:
%   DEM:        Struct containing information of the digital elevation
%               model at full 90m resolution, with fields:
%               .elevation_m:   Elevation data
%               .lat:           Latitude
%               .lon:           Longitude
%               .centroid_ID:   Only if centroids provided as input or if
%                               centroids input set to 1.
%   centroids:  Climada centroids struct with fields:
%               .elevation_m:   Elevation data
%               .lat:           Latitude
%               .lon:           Longitude
%               .centroid_ID:   Only if centroids provided as input or if
%                               centroids input set to 1.
%               .onLand:        Set to 0 if .elevation_m <0, 1 otherwise
%               .admin0_name    Country name
%               .admin0_ISO3    ISO 3 country code
% MODIFICATION HISTORY:
%   Gilles Stassen 20150107
%   Gilles Stassen 20150224     fixed some bugs in the plotting routines
%                               and added messages to warn of DEM edges
%   Gilles Stassen 20150225     cleanup and added automatic download and
%                               unzip feature
%-

DEM =[];

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if exist(climada_global.map_border_file, 'file')
    load(climada_global.map_border_file)
end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% gui to select srtm data if not provided
if ~exist('srtm_dir','var') || isempty(srtm_dir)
    srtm_dir = uigetdir(module_data_dir,'Choose srtm tile(s)');
end

if ~strcmp(srtm_dir, 'DL') && ~exist(srtm_dir, 'dir')
    cprintf([1 0 0], 'ERROR: Directory not found')
    return;
end

if ~exist('centroidsORcountry', 'var'),     centroidsORcountry  = [];   end
if ~exist('DEM_save_file',      'var'),     DEM_save_file       = [];   end
if ~exist('smooth',             'var'),     smooth              = [];   end
if ~exist('check_plot',         'var'),     check_plot          = 0;    end

if ~isempty(centroidsORcountry)
    if isstruct(centroidsORcountry)
        centroids   = centroidsORcountry; clear centroidsORcountry
        if isfield(centroids,'countryname')
            country_name    =   centroids.countryname;
            [country_name, ~, occurrence]=unique(country_name);
            country_name    =   country_name(mode(occurrence));
            country_name    =   country_name{1};
        else
            country_name = [];
        end
        rect        = [min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)];
    elseif numel(centroidsORcountry) == 4
        centroids   = [];
        rect        = centroidsORcountry; clear centroidsORcountry
    elseif ischar(centroidsORcountry)
        centroids   = [];
        country_name= centroidsORcountry; clear centroidsORcountry
        [country_name,country_ISO3,shape_index] = climada_country_name(country_name);
        % bb          = shapes(shape_index).BoundingBox;    % countries with colonies pose problems here...
        bb          = [min(shapes(shape_index).X) min(shapes(shape_index).Y)
                       max(shapes(shape_index).X) max(shapes(shape_index).Y)];
        rect        = [bb(:,1)' bb(:,2)']; clear bb
    end
else
    centroids   = []; clear centroidsORcountry
    [country_name,country_ISO3,shape_index] = climada_country_name('Single');
    if isempty(country_name), return; end % error message already printed in climada_country_name
    % bb          = shapes(shape_index).BoundingBox;    % countries with colonies pose problems here...
    bb          = [min(shapes(shape_index).X) min(shapes(shape_index).Y)
               max(shapes(shape_index).X) max(shapes(shape_index).Y)];
    rect        = [bb(:,1)' bb(:,2)']; clear bb
end
if ~isempty(country_name), cntry_str = sprintf(' for %s',country_name'); else cntry_str = ''; end

% load srtm tile from internet and unzip
if strcmp(srtm_dir, 'DL')
    dl_check = 1;
    clear srtm_dir
    
    % conversion to srtm tile indices
    srtm_min_lon_ndx = ceil(72 * (rect(1) + 180)/(179.28+180.00));
    srtm_max_lon_ndx = ceil(72 * (rect(2) + 180)/(179.28+180.00));
    srtm_min_lat_ndx = ceil(24 * (60 - rect(3)) /( 60.00+ 57.83));
    srtm_max_lat_ndx = ceil(24 * (60 - rect(4)) /( 60.00+ 57.83));
    
    % construct filenames
    n_tiles = (1+srtm_max_lon_ndx-srtm_min_lon_ndx)*(1+srtm_min_lat_ndx-srtm_max_lat_ndx);
    
    if n_tiles > 9
        warn_msg = sprintf('WARNING: Your specified region of interest requires %i DEM tiles. \n\t \t Computation may be slow and Matlab may crash. Are you sure you wish to continue? (y/n) ',n_tiles);
        response = input(warn_msg,'s');
        if ~strcmp(response,'y')
            fprintf('aborting\n')
            return
        end
    end   
    
    [I,J]   = meshgrid([srtm_min_lon_ndx: srtm_max_lon_ndx],[srtm_max_lat_ndx: srtm_min_lat_ndx]);
    t0 = clock;
    format_str = '%s';
    for tile_i = 1 : n_tiles
        srtm_fN {tile_i}    = strcat('srtm_',num2str(I(tile_i),'%02.0f'),'_',num2str(J(tile_i),'%02.0f'));
        srtm_dir{tile_i}    = [module_data_dir filesep 'system' filesep srtm_fN{tile_i}];
        srtm_URL{tile_i}    = ['ftp://srtm.csi.cgiar.org/SRTM_V41/SRTM_Data_GeoTiff/' srtm_fN{tile_i} '.zip'];
        
        if exist([srtm_dir{tile_i} filesep srtm_fN{tile_i} '.tif'],'file') && ...
                exist([srtm_dir{tile_i} filesep srtm_fN{tile_i} '.hdr'],'file')
            substr = sprintf('%s already exists - skipping', srtm_fN{tile_i});
            skip_file = 1;
        else
            % delete existing folder to avoid any unzipping issues
            if exist(srtm_dir{tile_i},'dir'), rmdir(srtm_dir{tile_i},'s'); end
            mkdir(srtm_dir{tile_i},'s');
            substr = sprintf('downloading and unzipping %s', srtm_fN{tile_i});
            skip_file = 0;
        end
        
        % progress management
        t_elapsed_tile   = etime(clock,t0)/tile_i;
        tiles_remaining  = n_tiles-tile_i;
        t_projected_sec   = t_elapsed_tile*tiles_remaining;
        if t_projected_sec<60
            msgstr = sprintf('%s, est. %3.0f sec left (%i/%i files)',substr, t_projected_sec, tile_i,n_tiles);
        else
            msgstr = sprintf('%s, est. %3.1f min left (%i/%i files)',substr, t_projected_sec/60,tile_i,n_tiles);
        end
        fprintf(format_str,msgstr); % write progress to stdout
        format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
        
        if ~skip_file
            unzip(srtm_URL{tile_i}, srtm_dir{tile_i})
        else
            pause(1)
        end
    end
else
    n_tiles = 1; I = 1; J = 1;
    tmp         = srtm_dir;     clear srtm_dir
    srtm_dir{1} = tmp;          clear tmp
end

fprintf('\nreading and processing DEM%s... ',cntry_str)

% struct containing the lat/lon extremes of each tile (in any order)
extremes.lon = [];
extremes.lat = [];

for tile_i = 1 : n_tiles
    srtm_files = dir(srtm_dir{tile_i});
    
    for file_i = 1 : numel(srtm_files)
        [~, ~, fE] = fileparts(srtm_files(file_i).name);
        
        if strcmp(fE, '.hdr')
            fid = fopen([srtm_dir{tile_i} filesep srtm_files(file_i).name]);
            scale_check = 0;
            while ~feof(fid),
                line = fgetl(fid);
                
                if scale_check
                    scale = str2num(line);
                    dlon = scale(1); dlat = scale(2);
                    scale_check =0;
                end
                if strfind(line,'ModelPixelScaleTag')
                    scale_check = 1;
                end
                
                if strfind(line,'Upper Left')
                    loc_i = strfind(line, '(');
                    loc_f = strfind(line, ')');
                    UL = str2num(line(loc_i+1:loc_f-1));
                end
                if strfind(line,'Lower Left')
                    loc_i = strfind(line, '(');
                    loc_f = strfind(line, ')');
                    LL = str2num(line(loc_i+1:loc_f-1));
                end
                if strfind(line,'Upper Right')
                    loc_i = strfind(line, '(');
                    loc_f = strfind(line, ')');
                    UR = str2num(line(loc_i+1:loc_f-1));
                end
                if strfind(line,'Lower Right')
                    loc_i = strfind(line, '(');
                    loc_f = strfind(line, ')');
                    LR = str2num(line(loc_i+1:loc_f-1));
                end
            end
            fclose(fid);
            extremes.lon = [extremes.lon UL(1) UR(1)];
            extremes.lat = [extremes.lat UL(2) LL(2)];
            break;
        end
    end
    
    for file_i = 1 : numel(srtm_files)
        [~, ~, fE] = fileparts(srtm_files(file_i).name);
        
        if strcmp(fE, '.tif')
            raw(I(tile_i),J(tile_i)).grid = imread([srtm_dir{tile_i} filesep srtm_files(file_i).name]);
            break;
        end
    end
end

DEM_grid = []; 

% Concatenate tiles
for i = srtm_min_lon_ndx: srtm_max_lon_ndx
    DEM_grid_j = [];
    for j = srtm_max_lat_ndx: srtm_min_lat_ndx
        DEM_grid_j = [DEM_grid_j ; raw(i,j).grid];
    end
    DEM_grid = [DEM_grid DEM_grid_j];
end
clear raw DEM_grid_j

reference_box = [min(extremes.lon) max(extremes.lon) min(extremes.lat) max(extremes.lat)];

DEM_grid(DEM_grid== min(min(DEM_grid))) = min(DEM_grid(DEM_grid~= min(min(DEM_grid))));
DEM_grid = double(DEM_grid);

% smooth the DEM if desired
if ~isempty(smooth) && any(smooth) && ~isnan(smooth)
    if isscalar(smooth)
        smooth_matrix = (1/smooth^2) .* ones(smooth);
    elseif ismatrix(smooth)
        smooth_matrix = smooth;
    end
    DEM_grid = filter2(smooth_matrix,DEM_grid);
end

% store as singleton arrays in DEM structure
[elev, lon, lat] = climada_grid2array(DEM_grid, reference_box);
fprintf('done \n')

if isstruct(centroids)
    % crop to rect
    lon_crop_ndx    = min(centroids.lon)<= lon & lon <= max(centroids.lon);
    lat_crop_ndx    = min(centroids.lat)<= lat & lat <= max(centroids.lat);
    lon_crop        = lon(lon_crop_ndx & lat_crop_ndx);
    lat_crop        = lat(lon_crop_ndx & lat_crop_ndx);
    elev_crop       = elev(lon_crop_ndx & lat_crop_ndx);
    
    DEM.elevation_m     = elev_crop';
    DEM.lon             = lon_crop';
    DEM.lat             = lat_crop';
    
    n_centroids = numel(centroids.centroid_ID);
    fprintf('processing centroid elevation... ');
    t0 = clock;
%     format_str = '%s';
%     for centroid_i = 1: n_centroids
%         
%         %r_i = climada_geo_distance(centroids.lon(centroid_i),centroids.lat(centroid_i),lon_crop,lat_crop);
%         r_i = sqrt((centroids.lon(centroid_i)-lon_crop).^2 + (centroids.lat(centroid_i)-lat_crop).^2);
%         [~,ndx] = min(r_i);
%         DEM.centroid_ID(ndx) = centroid_i;
%         centroids.elevation_m(centroid_i) = elev_crop(ndx)';
%         % the progress management
%         mod_step = 100;
%         if mod(centroid_i,mod_step)==0
%             t_elapsed_event   = etime(clock,t0)/centroid_i;
%             events_remaining  = n_centroids-centroid_i;
%             t_projected_sec   = t_elapsed_event*events_remaining;
%             if t_projected_sec<60
%                 msgstr = sprintf('est. %3.0f sec left (%i/%i centroids)',t_projected_sec,   centroid_i,n_centroids);
%             else
%                 msgstr = sprintf('est. %3.1f min left (%i/%i centroids)',t_projected_sec/60,centroid_i,n_centroids);
%             end
%             fprintf(format_str,msgstr); % write progress to stdout
%             format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
%         end
%     end

%     F_DEM = scatteredInterpolant(DEM.lon',DEM.lat',DEM.elevation_m');
    F_DEM = scatteredInterpolant(lon_crop,lat_crop,elev_crop);
    centroids.elevation_m = F_DEM(centroids.lon',centroids.lat')';
    fprintf('done \n')
%     fprintf(format_str,sprintf('processing centroid elevation%s took %3.0f seconds \n',cntry_str, etime(clock,t0)));
    
    if ~dl_check % only necessary when existing data is selected, avoid warning message due to imperfect cropping
        % deal with DEM edges if they do not reach extent of centroids
        if min(DEM.lon) > min(centroids.lon) || max(DEM.lon) < max(centroids.lon)
            cprintf([1 0.25 0.25], ['WARNING: longitudinal extent of centroids exceeds '...
                'that of DEM \n \t \t Elevation not processed for some centroids \n'])
        end
        if min(DEM.lat) > min(centroids.lat) || max(DEM.lat) < max(centroids.lat)
            cprintf([1 0.25 0.25], ['WARNING: latitudinal extent of centroids exceeds '...
                'that of DEM \n \t \t Elevation not processed for some centroids \n'])
        end
    end
    
    centroids.elevation_m(min(DEM.lon) > centroids.lon) = NaN;
    centroids.elevation_m(max(DEM.lon) < centroids.lon) = NaN;
    centroids.elevation_m(min(DEM.lat) > centroids.lat) = NaN;
    centroids.elevation_m(max(DEM.lat) < centroids.lat) = NaN;

    if ~exist('rect','var') || isempty(rect)
        rect = [min(DEM.lon) max(DEM.lon) min(DEM.lat) max(DEM.lat)];
    end
    
elseif ~isempty(rect)
        
    if min(DEM.lon) > rect(1) || max(DEM.lon) < rect(2)
        cprintf([1 0.25 0.25], ['WARNING: DEM does not cover longitudinal extent'...
            'defined by centroids rect \n \t \t Spatial extent of centroids limited to DEM \n'])
    end

    if min(DEM.lat) > rect(3) || max(DEM.lat) < rect(4)
        cprintf([1 0.25 0.25], ['WARNING: DEM does not cover latitudinal extent'...
            'defined by centroids rect \n \t \t Spatial extent of centroids limited to DEM \n'])
    end

    fprintf('generating centroids from DEM...');
    
    % crop to rect
    lon_crop_ndx    = rect(1) <= lon & lon <= rect(2);
    lat_crop_ndx    = rect(3) <= lat & lat <= rect(4);
    lon             = lon(lon_crop_ndx & lat_crop_ndx);
    lat             = lat(lon_crop_ndx & lat_crop_ndx);
    elev            = elev((lon_crop_ndx & lat_crop_ndx));
    
    DEM.elevation_m     = elev';
    DEM.lon             = lon';
    DEM.lat             = lat';
    
    % Generate centroids struct at same resolution as DEM if not provided
    centroids.lon           = lon';
    centroids.lat           = lat';
    centroids.elevation_m   = elev';
    n_centroids             = numel(centroids.lon);
    centroids.centroid_ID   = [1:n_centroids];
    DEM.centroid_ID         = [1:n_centroids];
    centroids.onLand        = ones(1,n_centroids);
    centroids.onLand(elev<0)= 0; % May be inaccurate when there are land points below sea level, but much faster than using inpolygon
    
    if exist('shapes','var')
        n_centroids = numel(centroids.lon);
        % accomodate for both climada global shp files, as well as shp
        % files downloaded from http://www.diva-gis.org/gdata
        if isfield(shapes,'NAME')  && length(shapes) == 1
            for i = 1 : n_centroids
                centroids.countryname{i} = shapes.NAME;
            end
            centroids.admin0_name = shapes.NAME;
        elseif exist('country_name','var')
            for i = 1 : n_centroids
                centroids.countryname{i} = country_name;
            end
            centroids.admin0_name = country_name;
        end
        if isfield(shapes,'ADM0_A3')  && length(shapes) == 1
            centroids.admin0_ISO3 = shapes.ADM0_A3;
        elseif exist('country_name','var')
            for i = 1 : n_centroids
                centroids.countryname{i} = country_ISO3;
            end
            centroids.admin0_name = country_ISO3;
        end
    end
    fprintf(' done \n');
else
    DEM.elevation_m = elev';
    DEM.lon         = lon';
    DEM.lat         = lat';
    rect            = reference_box;
end
clear elev lon lat

if ~isempty(DEM_save_file)
    save(DEM_save_file,'DEM', 'centroids');
end

if check_plot
    % relief plot
    figure('Name', '2D Relief Plot', 'color', 'w');
    hold on
    title(sprintf('Digital Elevation Model %s', cntry_str))

    if isempty(which('climada_DEM_plot'))
        fprintf('Download DEM plotting function from http://ch.mathworks.com/matlabcentral/fileexchange/36380-dem--shaded-relief-image-plot--digital-elevation-model-\n')
        fprintf('and rename dem -> climada_DEM_plot for awesome relief plot \n')
        fprintf('using imagesc instead\n')
        [s_x,s_y]   = size(DEM_grid);
        tmp_x       = linspace(reference_box(1),reference_box(2),s_x);
        tmp_y       = linspace(reference_box(3),reference_box(4),s_y);
        imagesc(tmp_x,tmp_y,DEM_grid)
    else
        [s_y,s_x]   = size(DEM_grid);
        tmp_x       = linspace(reference_box(1),reference_box(2),s_x);
        tmp_y       = linspace(reference_box(3),reference_box(4),s_y);
        climada_DEM_plot(tmp_x,fliplr(tmp_y),DEM_grid)
    end
    xlabel('Longitude');
    ylabel('Latitude');
    axis equal
    axis(rect)
    set(gca,'Ydir','normal')
    hold off
end

return

