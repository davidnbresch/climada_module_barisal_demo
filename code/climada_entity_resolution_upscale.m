function [entity_hr] = climada_entity_resolution_upscale(entity,x_factor,y_factor)
% climada upscale the resolution of a given entity struct
% NAME:
%   climada_entity_resolution_upscale
% PURPOSE:
%   given a lower resolution entity struct, this function interpolates the
%   existing entity assets to a new resolution defined by scaling factors
%   along the x and y directions. The new entity asset values are scaled
%   accordingly, ensuring that the total asset value is correct.
% CALLING SEQUENCE:
%   [entity_hr] = climada_entity_resolution_upscale(entity,x_factor,y_factor)
% EXAMPLE:
%   [entity_hr] = climada_entity_resolution_upscale(entity)
% INPUTS:
%   entity          : lower resolution entity structure, with entity.assets field
%   x_factor        : scaling factor for the x direction
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   entity_hr       : high resolution entity struct 
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com 20141104
%-

entity_hr = entity;

err_str = 'WARNING: One or more inputs not provided - entity not scaled';

if ~exist('entity','var'),      fprintf(err_str);   return;     end
if ~exist('x_factor','var'),    fprintf(err_str);   return;     end
if ~exist('y_factor','var'),    fprintf(err_str);   return;     end
if x_factor ==1 && y_factor==1,                     return;     end % routine unnecessary if scaling set to 1

% Store total asset value
total_value = sum(entity.assets.Value);
if isfield(entity.assets,'Value_today'),total_value_today = sum(entity.assets.Value_today);end

% Upscale lat lon resolution
[lon_hr, lat_hr] = resolution_upscale(entity.assets.Longitude,...
    entity.assets.Latitude,x_factor,y_factor);

entity_hr.assets.Longitude  = lon_hr;
entity_hr.assets.Latitude   = lat_hr;

% Interpolate values in each field accordingly
flds_i = fieldnames(entity.assets);
for i = 1 : numel(flds_i)
    fld = entity.assets.(flds_i{i});
    if numel(fld) ==numel(entity.assets.Longitude) && isnumeric(fld) &&...
            ~(strcmp(flds_i{i}, 'Longitude') || strcmp(flds_i{i}, 'Latitude'))
        entity_hr.assets.(flds_i{i}) = griddata(entity.assets.Longitude,...
            entity.assets.Latitude,entity.assets.(flds_i{i}),lon_hr,lat_hr);
        % Set griddata artifacts to zero
        entity_hr.assets.(flds_i{i})(entity_hr.assets.(flds_i{i}) < 1 |...
            isnan(entity_hr.assets.(flds_i{i})))=0;
    end
end



% Scale values appropriately
upscaled_value = sum(entity_hr.assets.Value);
entity_hr.assets.Value = entity_hr.assets.Value .*...
    (total_value / upscaled_value);

if isfield(entity.assets,'Value_today')
    entity_hr.assets.Value_today(entity_hr.assets.Value_today < 1 |...
        isnan(entity_hr.assets.Value_today))=0;
    upscaled_value_today = sum(entity_hr.assets.Value_today);
    entity_hr.assets.Value_today = entity_hr.assets.Value_today .*...
        (total_value_today / upscaled_value_today);
end
end