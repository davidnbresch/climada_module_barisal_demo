function [centroids, entity] = climada_clip_centroids_entity(...
    centroids, entity, bounding_box, centroids_scale_factor, entity_scale_factor, check_country)
% climada
% NAME:
%   climada_clip_centroids_entity
% PURPOSE:
%   Given low resolution centroids and entity assets on country level
%   (generated using climada_create_GDP_entity) this function uses the 2D  
%   linear interpolation methods employed in Matlab's griddata to produce a
%   high resolution (defined by the scale factors) set of centroids and
%   assets. As a result of the increase in asset resolution, and use of a 
%   high resolution border shape file, some assets will be located
%   offshore. This is an issue for storm surge, which is solved using the
%   climada_dist_nearest_coast function.
%
% CALLING SEQUENCE:
%   [centroids_HR, entity_HR] = climada_clip_centroids_entity(centroids_LR, entity_LR, bounding_box, [centroids_scale_factor], [entity_scale_factor])
% EXAMPLE:
%   [centroids_HR, entity_HR] = climada_clip_centroids_entity(centroids_LR, entity_LR, bounding_box)
% INPUTS:
%   centroids_LR:   The low resolution centroids created by climada_create_GDP_entity
%   entity_LR:      The low resolution entity created by climada_create_GDP_entity
%   bounding_box:   An array of size 4 bounding the region of interest,
%   defined by [min_lon max_lon min_lat max_lat]
% OPTIONAL INPUT PARAMETERS:
%   centroids_scale_factor:     Can be a single number (scales by same
%   factor along lat and lon), or an array of size 2, such that the first
%   index indicates the lon, and the second the lat scale factors. Default
%   value is set to 2 along both lat and lon.
%   entity_scale_factor:        Same concept as the centroids_scale_factor.
%   check_country:              Specify the country of interest. Since the
%   function was created for the Barisal module, default is set to
%   Bangladesh.
% OUTPUTS:
%   centroids_HR:     High resolution centroids
%   entity_HR:        High resolution entity assets
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com, 20141121
%   Gilles Stassen, gillesstassen@hotmail.com, 20141218 change variable
%                       names whole_world_borders.lon/lat -> shapes.X/Y
%   Gilles Stassen, gillesstassen@hotmail.com, 20141223 add check_country input arg.
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('centroids', 'var'),              centroids               = [];     end
if ~exist('entity', 'var'),                 entity                  = [];     end
if ~exist('bounding_box', 'var'),           bounding_box            = [];     end
if ~exist('centroids_scale_factor', 'var'), centroids_scale_factor  = [];     end
if ~exist('entity_scale_factor', 'var'),    entity_scale_factor     = [];     end
if ~exist('check_country','var'),           check_country           = 'Bangladesh'; end

if isempty(centroids) && isempty(entity)
    fprintf('ERROR: Please provide an input of either centroids, entity or both');
    return;
end

% If scale factor not specified, use default resolution increase 2x2
if isempty(centroids_scale_factor),c_s_f_x = 2; c_s_f_y = 2;end
if isempty(entity_scale_factor),e_s_f_x = 2; e_s_f_y = 2;end

% If scale_factor is an array of two (x & y) values, scale x & y coords
% independently, otherwise, use same factor for both.
if ~isempty(centroids_scale_factor)
    switch numel(centroids_scale_factor)
        case 2
            c_s_f_x = centroids_scale_factor(1);
            c_s_f_y = centroids_scale_factor(2);
        case 1
            c_s_f_x = centroids_scale_factor;
            c_s_f_y = centroids_scale_factor;
    end
end
if ~isempty(entity_scale_factor)
    switch numel(entity_scale_factor)
        case 2
            e_s_f_x = entity_scale_factor(1);
            e_s_f_y = entity_scale_factor(2);
        case 1
            e_s_f_x = entity_scale_factor;
            e_s_f_y = entity_scale_factor;
    end
end

if numel(bounding_box) ~= 4
    fprintf('ERROR: Please provide a bounding box defining the clipping area');
    return;
end
% Set to 1 to move waterborne entity assets onto coast. Uses inpolygon, so
% can be slow for high res borders. For the damage calculation, we take the
% only centroid geolocations are used, and actual asset lat lon is
% irrelevant, other than in the process of encoding assets to centroids...
move_assets_to_land = 1; 

load(climada_global.map_border_file);
[~,~,shape_index] = climada_country_name(check_country);
shapes = shapes(shape_index);

if ~isempty(centroids)
    % Trim centroids struct
    tmp_lon = centroids.lon < bounding_box(2) & ...
        centroids.lon > bounding_box(1);
    tmp_lat = centroids.lat< bounding_box(4) & ...
        centroids.lat > bounding_box(3);
    
    % For each field with the same number of data points as there are
    % centroids, trim to bounding box
    flds = fieldnames(centroids);
    no_c_ori = numel(centroids.centroid_ID);
    for i = 1 : numel(flds)
        if numel(centroids.(flds{i})) == no_c_ori
            centroids.(flds{i}) = centroids.(flds{i})((tmp_lon & tmp_lat));
        end
    end
    
    % Increase resolution of centroids by linear interpolation (using griddata)
    centroids = climada_centroids_resolution_upscale(centroids,c_s_f_x,c_s_f_y);
    
    % Reassign logical onLand values based on higher resolution Bangladesh
    % polygon and centroids
    in = inpolygon(centroids.lon,centroids.lat,shapes.X,shapes.Y);
    centroids.onLand = false(size(centroids.centroid_ID));
    centroids.onLand(in) = 1;
end

if ~isempty(entity)
    
    tmp_lon = entity.assets.lon < bounding_box(2) & ...
        entity.assets.lon > bounding_box(1);
    tmp_lat = entity.assets.lat < bounding_box(4) & ...
        entity.assets.lat > bounding_box(3);
    
    % Trim entity structure to bounding box
    flds = fieldnames(entity.assets);
    no_ea_ori = numel(entity.assets.Value);
    for i = 1 : numel(flds)
        if numel(entity.assets.(flds{i})) == no_ea_ori
            entity.assets.(flds{i}) = entity.assets.(flds{i})((tmp_lon & tmp_lat));
        end
    end
    
    entity = climada_entity_resolution_upscale(entity,e_s_f_x,e_s_f_y);

    % Move newly created asset points from water to land, if they are
    % waterborne w.r.t. new, high res borders, and as a result of upscaling
    if move_assets_to_land
        
        fprintf('moving assets onto land...')
        
        inval = entity.assets.Value > 0;
        tmp_lon = entity.assets.lon(inval);
        tmp_lat = entity.assets.lat(inval);
        
        inland = inpolygon(tmp_lon,tmp_lat,shapes.X,shapes.Y); % can take some time
        
        tmp_lon_water = tmp_lon(~inland);
        tmp_lat_water = tmp_lat(~inland);

        [tmp_lon_water, tmp_lat_water] = climada_dist_nearest_coast(tmp_lon_water,tmp_lat_water);

        tmp_lon(~inland) = tmp_lon_water;
        tmp_lat(~inland) = tmp_lat_water;
        
        entity.assets.lon(inval) = tmp_lon;
        entity.assets.lat(inval) = tmp_lat;
        
        fprintf(' done \n')
    end
end
