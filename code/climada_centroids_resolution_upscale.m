function [centroids_hr] = climada_centroids_resolution_upscale(centroids,x_factor,y_factor)

% climada upscale the resolution of a given centroids struct
% NAME:
%   climada_centroids_resolution_upscale
% PURPOSE:
%   given a lower resolution centroids struct, this function interpolates the
%   existing centroids to a new resolution defined by scaling factors
%   along the x and y directions.
% CALLING SEQUENCE:
%   [centroids_hr] = climada_centroids_resolution_upscale(centroids,x_factor,y_factor)
% EXAMPLE:
%   [centroids_hr] = climada_centroids_resolution_upscale(centroids)
% INPUTS:
%   centroids       : lower resolution entity structure, with entity.assets field
%   x_factor        : scaling factor for the x direction
%   y_factor        : scaling factor for the y direction
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   centroids_hr       : high resolution centroids struct
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com 20141104
%-

centroids_hr = centroids;

err_str = 'WARNING: One or more inputs not provided - centroids not scaled';

if ~exist('centroids','var'),   fprintf(err_str);   return;     end
if ~exist('x_factor','var'),    fprintf(err_str);   return;     end
if ~exist('y_factor','var'),    fprintf(err_str);   return;     end

% Scale centroids lat lon
[lon_hr, lat_hr] = resolution_upscale(centroids.Longitude,...
    centroids.Latitude,x_factor,y_factor);

centroids_hr.Longitude = lon_hr;
centroids_hr.Latitude = lat_hr;

% Interpolate the data in each field of the centroids struct
flds_i = fieldnames(centroids);
for i = 1 : numel(flds_i)
    
    fld = centroids.(flds_i{i});
    if numel(fld) > 1 && isnumeric(fld) &&...
            ~(strcmp(flds_i{i}, 'Longitude') || strcmp(flds_i{i}, 'Latitude'))
        centroids_hr.(flds_i{i}) = griddata(centroids.Longitude,...
            centroids.Latitude,centroids.(flds_i{i}),lon_hr,lat_hr);
    end
end

% Reassign centroid_IDs
centroids_hr.centroid_ID = linspace(1,numel(centroids_hr.centroid_ID),numel(centroids_hr.centroid_ID));

end