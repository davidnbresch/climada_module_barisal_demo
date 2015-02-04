function [c_lon, c_lat, d_lon, d_lat] = climada_dist_nearest_coast(lon, lat, check_country, polygon)
% climada
% NAME:
%   climada_dist_nearest_coast
% PURPOSE:
%   This function is used in climada_clip_centroids_entity to move
%   waterborne entity assets onto the nearest shore. It may be used in a
%   more general setting, moving a given point onto the nearest edge of a
%   polygon.
% CALLING SEQUENCE:
%   [c_lon, c_lat, d_lon, d_lat] = climada_dist_nearest_coast(lon, lat, [polygon])
% EXAMPLE:
%   [c_lon, c_lat, d_lon, d_lat] = climada_dist_nearest_coast(lon, lat)
% INPUTS:
%   lon:   The longitude of a single point
%   lat:   The latitude of a single point
% OPTIONAL INPUT PARAMETERS:
%   polygon:    A polygon structure with .lat and .lon fields. If not
%               provided, the function uses the border polygon from
%               climada_global.map_border_file.
% OUTPUTS:
%   c_lon, c_lat:     Lon & lat of nearest point on coast (polygon)
%   d_lon, d_lat:     Difference in lon & lat from coast to input point.
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com, 20141201
%   Gilles Stassen, gillesstassen@hotmail.com, 20141208 accept vector for lon & lat input
%   Gilles Stassen, gillesstassen@hotmail.com, 20141218 change variable
%                       names whole_world_borders.lon/lat -> shapes.X/Y
%   Gilles Stassen, gillesstassen@hotmail.com, 20141223 add check_country
%   Gilles Stassen, 20150112 incorporate climada_geo_distance
%-

% Init
c_lon = []; c_lat = []; d_lon = []; d_lat = [];

% Check arguments
if ~exist('polygon', 'var'),    polygon = [];                   end
if ~exist('check_country','var'), check_country = 'Bangladesh'; end % Since created for Barisal module
if ~exist('lon', 'var') || ~exist('lat', 'var'),    return;     end
if isempty(lon) || isempty(lat),                    return;     end
if numel(lon) ~= numel(lat)
    fprintf('ERROR: Size of lat and lon must match. Unable to proceed.')
    return;
end

global climada_global

if isfield(climada_global,'map_border_file') && isempty(polygon)
    load(climada_global.map_border_file)
    [~,~,shape_index] = climada_country_name(check_country);
    polygon.lon = shapes(shape_index).X;
    polygon.lat = shapes(shape_index).Y;
elseif ~exist(shapes,'var') && isempty(polygon)
    polygon = input('ERROR: No coastal information available, please provide the variable name of a polygon structure, with fields "lat" and "lon":');
end

% Does this work???
for ndx = 1 : numel(lon)
    % r = sqrt((polygon.lon - lon(ndx)).^2 + (polygon.lat - lat(ndx)).^2);
    r = climada_geo_distance(lon(ndx),lat(ndx),polygon.lon,polygon.lat);
    if sum(r == min(r))==1
        c_lon(ndx) = polygon.lon(r == min(r));
        c_lat(ndx) = polygon.lat(r == min(r));
    else
        tmp = polygon.lon(r == min(r));
        c_lon(ndx) = tmp(1);
        tmp = polygon.lat(r == min(r));
        c_lat(ndx) = tmp(1);
    end
end
d_lon = c_lon - lon;
d_lat = c_lat - lat;
return; % TEST


% Determine resolution of polygon vertices (assumes they are uniformly
% distributed. If not, this is a good enough approximation for current
% purposes)
res_poly_lon = (max(unique(polygon.lon)) - min(unique(polygon.lon)))/numel(unique(polygon.lon)); % Longitude
res_poly_lat = (max(unique(polygon.lat)) - min(unique(polygon.lat)))/numel(unique(polygon.lat)); % Latitude

poly_crop.lon = [];
poly_crop.lat = [];

% Create box of width i * res_poly_lon and length i * res_poly_lat,
% centered on input point. Box size increases on every iteration, until it
% encloses a portion of the polygon - the enclosed polygon is poly_crop.
for ndx=1: numel(lon)
    i = 1;
    while (isempty(poly_crop.lon) || isempty(poly_crop.lat))
        min_lon = lon(ndx) - i*res_poly_lon;
        max_lon = lon(ndx) + i*res_poly_lon;
        
        min_lat = lat(ndx) - i*res_poly_lat;
        max_lat = lat(ndx) + i*res_poly_lat;
        
        poly_crop.lon = polygon.lon(((min_lon <= polygon.lon) & (polygon.lon <= max_lon)) &...
            ((min_lat <= polygon.lat) & (polygon.lat <= max_lat)));
        poly_crop.lat = polygon.lat(((min_lon <= polygon.lon) & (polygon.lon <= max_lon)) &...
            ((min_lat <= polygon.lat) & (polygon.lat <= max_lat)));
        
        i = i + 1;
    end
    r = sqrt((poly_crop.lon - lon(ndx)).^2 + (poly_crop.lat - lat(ndx)).^2);
    
    if sum(r == min(r))==1
        c_lon(ndx) = poly_crop.lon(r == min(r));
        c_lat(ndx) = poly_crop.lat(r == min(r));
    else
        tmp = poly_crop.lon(r == min(r));
        c_lon(ndx) = tmp(1);
        tmp = poly_crop.lat(r == min(r));
        c_lat(ndx) = tmp(1);
    end
end
d_lon = c_lon - lon;
d_lat = c_lat - lat;

end