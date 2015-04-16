function centroids = climada_generate_centroids(centroids_rectORcountry_name, resolution_km, buffer_check, check_plot)
% climada generate high resolution centroids
% MODULE:
%   barisal_demo
% NAME:
%   climada_generate_HR_centroids
% PURPOSE:
%   Given a rectangle defining the location of interest, generate an evenly
%   spaced rectilinear grid of hazard centroids
% CALLING SEQUENCE:
%   centroids = climada_generate_HR_centroids(centroids_rect, resolution_km);
%   centroids = climada_generate_HR_centroids([min_lon max_lon min_lat max_lat], resolution_km);
% EXAMPLE:
%   centroids = climada_generate_HR_centroids([75.0 77.8 21.3 25.1], resolution_km);
% INPUTS:
%   centroids_rect: 4-element row vector defining the longitude and latitude
%                   limits of the study region
% OPTIONAL INPUT PARAMETERS:
%   resolution_km:  specify the centroid resolution (default = 1 km)
% OUTPUTS:
%   centroids:      climada centroids struct with fields
%                     .Longitude
%                     .Latitude
%                     .onLand
%                     .centroid_ID
%                     .countryname - cell array same size as .centroid_ID
%                     .admin0_name - country name char array
%                     .admin0_ISO3 - ISO 3 country code
%                     .comment
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150119
% Gilles Stassen, 20150128, added .comment field
% Gilles Stassen, 20150326, added buffer
% Gilles Stassen, 20150408, increased buffer size
%-
centroids = [];

global climada_global
if ~climada_init_vars, return; end

if exist(climada_global.map_border_file,'file')
    load(climada_global.map_border_file)
    shapes_check = 1;
else
    shapes_check = 0;
end

country_check = 0;
if ~exist('centroids_rectORcountry_name','var')
    if ~shapes_check
        % need shapes for country input
        fprintf('ERROR: no world border info found \n')
        return
    end
    
    if length(shapes) == 1 && isfield(shapes,'BoundingBox')
        centroids_rect =[shapes.BoundingBox(:,1)' shapes.BoundingBox(:,2)'];
    else
        [country_name,ISO_3,shape_index] = climada_country_name('Single');
        if isempty(country_name), return; end % error message already printed in climada_country_name
        bb             = [min(shapes(shape_index).X) min(shapes(shape_index).Y)
                          max(shapes(shape_index).X) max(shapes(shape_index).Y)];
        centroids_rect = [bb(:,1)' bb(:,2)']; clear bb
        country_check  = 1;
    end
    
elseif ischar(centroids_rectORcountry_name) 
    if ~shapes_check
        % need shapes for country input
        fprintf('ERROR: no world border info found \n')
        return
    end
    %input is country name
    country_name = centroids_rectORcountry_name; clear centroids_rectORcountry_name
    [country_name,ISO_3,shape_index] = climada_country_name(country_name);
    if isempty(country_name),   return;             end % error message already printed in climada_country_name
    if length(shapes) == 1,     shape_index = 1;    end
    bb             = [min(shapes(shape_index).X) min(shapes(shape_index).Y)
                      max(shapes(shape_index).X) max(shapes(shape_index).Y)];
    centroids_rect = [bb(:,1)' bb(:,2)']; clear bb
    country_check  = 1;
elseif isnumeric(centroids_rectORcountry_name) && length(centroids_rectORcountry_name) == 4
    % input is centroids rect
    centroids_rect = centroids_rectORcountry_name; clear centroids_rectORcountry_name
else
    fprintf('ERROR: invalid input argument \n')
    return
end

if ~exist('resolution_km','var') || isempty(resolution_km)
    fprintf('WARNING: no resolution specified, resorting to default 1 km \n')
    resolution_km = 1;
end
resolution_ang = resolution_km / (111.12);

if ~exist('buffer_check',   'var'),     buffer_check    = 1;    end
if ~exist('check_plot',     'var'),     check_plot      = 0;    end

min_lon = centroids_rect(1);
max_lon = centroids_rect(2);
min_lat = centroids_rect(3);
max_lat = centroids_rect(4);

n_lon = round((max_lon - min_lon)/(resolution_ang)) + 1;
n_lat = round((max_lat - min_lat)/(resolution_ang)) + 1;

fprintf(sprintf('generating centroids at %3.2f km resolution... ', resolution_km));
t0 = clock;
% construct regular grid
for i = 0 : n_lon - 1
    ndx = i * n_lat;
    centroids.lat(1,ndx + 1 : ndx + n_lat)= (1:n_lat) .* resolution_ang + min_lat;
    centroids.lon(1,ndx + 1 : ndx + n_lat)= (n_lon - i) .* resolution_ang + min_lon;
end

if shapes_check
    if buffer_check
        % take only high resolution centroids for country
        country_ndx   = inpolygon(centroids.lon,centroids.lat,shapes(shape_index).X,shapes(shape_index).Y);
        centroids.lon = centroids.lon(country_ndx);
        centroids.lat = centroids.lat(country_ndx);
        
        % generate coarser grid for buffer centroids (resolution reduced by
        % factor of 4)
        buffer_resolution_ang = resolution_ang * 4;
        n_lon = round((max_lon - min_lon)/(buffer_resolution_ang)) + 2;
        n_lat = round((max_lat - min_lat)/(buffer_resolution_ang)) + 3;

        for i = 0 : n_lon + 1
            ndx = i * n_lat;
            buffer.lat(1,ndx + 1 : ndx + n_lat)= (-1 : n_lat-2) .* buffer_resolution_ang + min_lat;
            buffer.lon(1,ndx + 1 : ndx + n_lat)= (n_lon - i) .* buffer_resolution_ang + min_lon;
        end
        
        % take only buffer centroids outside country border
        buffer_ndx = inpolygon(buffer.lon,buffer.lat,shapes(shape_index).X,shapes(shape_index).Y);
        buffer.lon = buffer.lon(~buffer_ndx);
        buffer.lat = buffer.lat(~buffer_ndx);
        
        % concatenate
        buffer_logical      = [zeros(size(centroids.lon)) ones(size(buffer.lon))];
        centroids.lon       = [centroids.lon buffer.lon];
        centroids.lat       = [centroids.lat buffer.lat];
    end
    
    
    if exist('country_name','var') 
        for i = 1 : length(centroids.lon)
            centroids.country_name{i} = country_name;
        end
        centroids.admin0_ISO3 = ISO_3;
    else
        if shapes_check && isfield(shapes,'NAME')
            for i = 1 : length(centroids.lon)
                centroids.country_name{i} = shapes.NAME;
            end
            centroids.admin0_name = shapes.NAME;
        end
        if isfield(shapes,'ADM0_A3')
            centroids.admin0_ISO3 = shapes.ADM0_A3;
        end
    end
end

centroids.centroid_ID = 1:length(centroids.lon);

% coastline is terrible.
centroids.onLand = ones(size(centroids.centroid_ID));
centroids.onLand(buffer_logical ==1) = 0;
% 
% coastline = climada_coastline_read;
% 
% if ~isempty(coastline)
%     in = inpolygon(centroids.lon,centroids.lat,coastline.lon,coastline.lat);
%     centroids.onLand        = true(size(centroids.centroid_ID));
%     centroids.onLand(~in & buffer_logical)   = false;
% else
%     fprintf('WARNING: coastline info not found, .onLand set to NaN \n')
% end
fprintf('done \n')
tf = clock;
fprintf('generating %i centroids for %s took %3.2f seconds \n',length(centroids.centroid_ID),country_name, etime(tf,t0));
centroids.comment = sprintf('%3.2f km resolution centroids, created on %s', resolution_km,datestr(now,'dd/mm/yyyy'));

if check_plot
    figure('name','Centroids','color','w')
    hold on
    climada_plot_world_borders
    scatter(centroids.lon(centroids.onLand ==1),centroids.lat(centroids.onLand == 1),'xr');
    scatter(centroids.lon(centroids.onLand ==0),centroids.lat(centroids.onLand == 0),'ob');
    if exist('country_name','var') && ~isempty(country_name)
        title(sprintf('Centroids for %s',country_name));
    end
    axis([min(centroids.lon) max(centroids.lon) min(centroids.lat) max(centroids.lat)])
    xlabel('Longitude')
    ylabel('Latitude')
end

return
