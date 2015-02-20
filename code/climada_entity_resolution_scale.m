function [entity] = climada_entity_resolution_scale(entity,x_factor,y_factor)
% climada upscale the resolution of a given entity struct
% NAME:
%   climada_entity_resolution_scale
% PURPOSE:
%   given an entity struct, this function interpolates the
%   existing entity assets to a new resolution defined by scaling factors
%   along the x and y directions. The new entity asset values are scaled
%   accordingly, ensuring that the total asset value is correct.
% CALLING SEQUENCE:
%   [entity] = climada_entity_resolution_upscale(entity,x_factor,y_factor)
% EXAMPLE:
%   [entity] = climada_entity_resolution_scale(entity)
% INPUTS:
%   entity          : entity structure, with entity.assets field
%   x_factor        : scaling factor for the x direction
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   entity          : scaled resolution entity struct 
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com 20141104
%   Gilles Stassen, gillesstassen@hotmail.com 20150220 renamed function ..._upscale -> ..._scale
%-

err_str = 'WARNING: One or more inputs not provided - entity not scaled';

if ~exist('entity','var'),      entity=climada_entity_load;     end
if ~exist('x_factor','var'),    fprintf(err_str);   return;     end
if ~exist('y_factor','var'),    fprintf(err_str);   return;     end
if x_factor ==1 && y_factor==1,                     return;     end % routine unnecessary if scaling set to 1

% Store total asset value
total_value = sum(entity.assets.Value);
if isfield(entity.assets,'Value_today'),total_value_today = sum(entity.assets.Value_today);end

% Scale lat lon resolution
[lon, lat] = resolution_scale(entity.assets.lon,...
    entity.assets.lat,x_factor,y_factor);

entity.assets.lon  = lon;
entity.assets.lat  = lat;

% Interpolate values in each field accordingly
flds_i = fieldnames(entity.assets);
for i = 1 : numel(flds_i)
    fld = entity.assets.(flds_i{i});
    if numel(fld) ==numel(entity.assets.lon) && isnumeric(fld) &&...
            ~(strcmp(flds_i{i}, 'lon') || strcmp(flds_i{i}, 'lat'))
        entity.assets.(flds_i{i}) = griddata(entity.assets.lon,...
            entity.assets.lat,entity.assets.(flds_i{i}),lon,lat);
        % Set griddata artifacts to zero
        entity.assets.(flds_i{i})(entity.assets.(flds_i{i}) < 1 |...
            isnan(entity.assets.(flds_i{i})))=0;
    end
end

% Scale values appropriately
scaled_value = sum(entity.assets.Value);
entity.assets.Value = entity.assets.Value .* (total_value / scaled_value);

if isfield(entity.assets,'Value_today')
    entity.assets.Value_today(entity.assets.Value_today < 1 |...
        isnan(entity.assets.Value_today))=0;
    scaled_value_today = sum(entity.assets.Value_today);
    entity.assets.Value_today = entity.assets.Value_today .*...
        (total_value_today / scaled_value_today);
end