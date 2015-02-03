function [DEM, centroids] = climada_read_srtm_DEM(srtm_dir, centroids, DEM_save_file, smooth, check_plot)
% climada
% MODULE:
%   barisal_demo
% NAME:
%   climada_read_srtm_DEM
% PURPOSE:
%   Read the digital elevation model data from the files in an existing
%   srtm directory. Data can be downloaded from http://srtm.csi.cgiar.org/SELECTION/inputCoord.asp
% CALLING SEQUENCE:
%   DEM = climada_read_srtm_DEM(srtm_dir, centroids, check_plot)
% EXAMPLE:
%   DEM = climada_read_srtm_DEM(srtm_dir, [], 1)
%   DEM = climada_read_srtm_DEM
% INPUTS:
%   srtm_dir:   The directory of an srtm data tile folder, containing at
%               least a .hdr and a .tif file
% OPTIONAL INPUT PARAMETERS:
%   centroids:  If centroids are provided as an input, the DEM will contain
%               elevation data sampled at the location of centroids, and
%               hence will have an extra field .centroid_ID. The extra
%               field is required for tc_surge_hazard_create if you wish to
%               provide your own topography data (e.g. srtm).
%               If this input is set to 1 (i.e. not a centroids struct), a
%               centroids struct will be generated at the same resolution
%               as the DEM.
%               If set to a 4-element vector (centroids_rect), these 4
%               points will define the area of interest, which is
%               subsequently cropped out.
%               NOTE:   It is only sensible to provide a centroids struct
%               as input if its resolution is significantly lower than that
%               of the DEM, otherwise, it is much faster to generate
%               centroids directly from the DEM.
%   smooth:     Can either be set to an integer N (smooth by default filter
%               specified by a matrix size NxN with values 1/N) or a
%               smoothing filter. Default = [] (no smoothing).
%   check plot: Specify whether to plot a relief of the DEM, default = 0
% OUTPUTS:
%   DEM:        Struct containing information of the digital elevation
%               model, with fields:
%               .elevation_m:   Elevation data
%               .lat:           Latitude
%               .lon:           Longitude
%               .centroid_ID:   Only if centroids provided as input or if
%                               centroids input set to 1.
% MODIFICATION HISTORY:
%   Gilles Stassen 20150107
%-

DEM =[];

global climada_global
if ~climada_init_vars,return;end % init/import global variables

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

if ~exist('srtm_dir','var')
    srtm_dir = uigetdir(module_data_dir,'Choose srtm tile');
end

if ~exist(srtm_dir, 'dir')
    fprintf('ERROR: Directory not found')
    return;
end

if ~exist('centroids', 'var'), centroids = []; end
if ~exist('DEM_save_file','var'), DEM_save_file = []; end
if ~exist('smooth','var'), smooth = []; end
if ~exist('check_plot', 'var'), check_plot = 0; end

srtm_files = dir(srtm_dir);

for file_i = 1 : numel(srtm_files)
    [~, ~, fE] = fileparts(srtm_files(file_i).name);
    
    if strcmp(fE, '.hdr')
        fid = fopen([srtm_dir filesep srtm_files(file_i).name]);
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
        reference_box = [UL(1) UR(1) LL(2) UL(2)];
        break;
    end
    
end

for file_i = 1 : numel(srtm_files)
    [~, ~, fE] = fileparts(srtm_files(file_i).name);
    
    if strcmp(fE, '.tif')
        DEM_grid = imread([srtm_dir filesep srtm_files(file_i).name]);
        DEM_grid(DEM_grid== min(min(DEM_grid))) = min(DEM_grid(DEM_grid~= min(min(DEM_grid))));
        DEM_grid = double(DEM_grid);
        
        if ~isempty(smooth) && any(smooth) && ~isnan(smooth)
            if isscalar(smooth)
                smooth_matrix = (1/smooth^2) .* ones(smooth);
            elseif ismatrix(smooth)
                smooth_matrix = smooth;
            end
            DEM_grid = filter2(smooth_matrix,DEM_grid);
        end
        [elev, lon, lat] = climada_grid2array(DEM_grid, reference_box);
        % elev(elev== min(elev)) = min(elev(elev~= min(elev)));
        break;
    end
end

if isstruct(centroids)
    
    lon_crop_ndx  = min(centroids.lon)<= lon & lon <= max(centroids.lon);
    lat_crop_ndx  = min(centroids.lat)<= lat & lat <= max(centroids.lat);
    lon_crop = lon(lon_crop_ndx & lat_crop_ndx);
    lat_crop = lat(lon_crop_ndx & lat_crop_ndx);
    elev_crop = elev(lon_crop_ndx & lat_crop_ndx);
    
    n_centroids = numel(centroids.centroid_ID);
    fprintf('processing centroid elevation \n');
    t0 = clock;
    format_str = '%s';
    for centroid_i = 1: n_centroids
        %r_i = climada_geo_distance(centroids.lon(centroid_i),centroids.lat(centroid_i),lon_crop,lat_crop);
        r_i = sqrt((centroids.lon(centroid_i)-lon_crop).^2 + (centroids.lat(centroid_i)-lat_crop).^2);
        [~,ndx] = min(r_i);
        DEM.elevation_m(centroid_i) = elev_crop(ndx);
        
        % the progress management
        mod_step = 100;
        if mod(centroid_i,mod_step)==0
            t_elapsed_event   = etime(clock,t0)/centroid_i;
            events_remaining  = n_centroids-centroid_i;
            t_projected_sec   = t_elapsed_event*events_remaining;
            if t_projected_sec<60
                msgstr = sprintf('est. %3.0f sec left (%i/%i centroids)',t_projected_sec,   centroid_i,n_centroids);
            else
                msgstr = sprintf('est. %3.1f min left (%i/%i centroids)',t_projected_sec/60,centroid_i,n_centroids);
            end
            fprintf(format_str,msgstr); % write progress to stdout
            format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
        end
    end
    DEM.centroid_ID = centroids.centroid_ID;
    DEM.lon         = centroids.lon;
    DEM.lat         = centroids.lat;
    fprintf(format_str,sprintf('processing DEM took %3.0f seconds \n',etime(clock,t0)));
    
    rect = [min(DEM.lon) max(DEM.lon) min(DEM.lat) max(DEM.lat)];
    
    %     DEM.elevation_m = griddata(lon_crop, lat_crop, elev_crop, centroids.lon, centroids.lat);
    %     DEM.centroid_ID = centroids.centroid_ID;
else
    if ~isempty(centroids)
        fprintf('generating centroids from DEM...');
        if numel(centroids) == 4 % centroids_rect has been provided as input
            rect = centroids;
            lon_crop_ndx  = rect(1) <= lon & lon <= rect(2);
            lat_crop_ndx  = rect(3) <= lat & lat <= rect(4);
            lon = lon(lon_crop_ndx & lat_crop_ndx);
            lat = lat(lon_crop_ndx & lat_crop_ndx);
            elev = elev((lon_crop_ndx & lat_crop_ndx));
        end
        clear centroids
        
        DEM.elevation_m = elev;
        DEM.lon         = lon;
        DEM.lat         = lat;
        
        % Generate centroids struct at same resolution as DEM if not provided
        centroids.lon     = lon;
        centroids.lat      = lat;
        centroids.elevation_m   = elev;
        n_centroids = numel(centroids.lon);
        centroids.centroid_ID   = [1:n_centroids]';
        DEM.centroid_ID         = [1:n_centroids]';
        centroids.onLand        = ones(n_centroids,1);
        centroids.onLand(elev<0)= 0; % May be inaccurate when there are land points below sea level, but much faster than using inpolygon
        
        if exist(climada_global.map_border_file,'file')
            load(climada_global.map_border_file)
            %     in = inpolygon(centroids.lon,centroids.lat,shapes.X,shapes.Y);
            %     centroids.onLand = false(size(centroids.centroid_ID));
            %     centroids.onLand(in) = 1;
            
            n_centroids = numel(centroids.lon);
            
            if isfield(shapes,'NAME')
                for i = 1 : n_centroids
                    centroids.countryname{i} = shapes.NAME;
                end
                centroids.admin0_name = shapes.NAME;
            end
            if isfield(shapes,'ADM0_A3')
                centroids.admin0_ISO3 = shapes.ADM0_A3;
            end
        end
        fprintf(' done \n');
    else
        DEM.elevation_m = elev;
        DEM.lon         = lon;
        DEM.lat         = lat;
        rect            = reference_box;
    end
end

if ~isempty(DEM_save_file)
    save(DEM_save_file,'DEM', 'centroids', 'DEM_grid');
end

if check_plot
    color_map = [
        40  54  154
        0   201  50
        30  211 104
        94  224 116
        162 235 130
        223 248 146
        246 229 149
        200 178 118
        162 126  94
        143  97  84
        162 125 116
        178 150 139
        199 176 170
        219 205 202
        236 228 226
        255 255 255
        ]./255;
    
    x = DEM.lon; y = DEM.lat;
    
    tmp_x = unique(x);
    tmp_y = unique(y);
    
    for i = 1 : numel(tmp_x)
        Z(i,:) = DEM.elevation_m(DEM.lon == tmp_x(i));
    end
    Z = Z';
    
    figure('Name', '2D Relief Plot', 'color', 'w');
    imagesc(tmp_x,tmp_y,Z);
    hold on
    axis(rect)
    axis equal
    set(gca,'YDir','reverse');
    colormap(color_map);
    caxis([-1 30]);
    xlabel('Longitude');
    ylabel('Latitude');
    hold off
    
    figure('Name', '3D Surface Plot', 'color', 'w');
    tri = delaunay(DEM.lon, DEM.lat);
    trisurf(tri,DEM.lon,DEM.lat,DEM.elevation_m);
    shading interp
    hold on
    axis(rect)
    set(gca,'YDir','reverse');
    colormap(color_map);
    caxis([-1 30]);
    xlabel('Longitude');
    ylabel('Latitude');
    view(330,70)
    set(gca,'ztick',[], 'zcolor', 'w');
    grid off 
    box off
    hold off
    
end

return

