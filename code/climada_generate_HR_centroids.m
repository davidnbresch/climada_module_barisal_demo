function centroids = climada_generate_HR_centroids(centroids_rect, resolution_km)

centroids = [];

global climada_global
if ~climada_init_vars, return; end

err_msg = 'ERROR: valid bounding box not specified';

if ~exist('centroids_rect','var')
    if exist(climada_global.map_border_file,'file')
        load(climada_global.map_border_file)
        if isfield(shapes,'BoundingBox')
            centroids_rect =[shapes.BoundingBox(:,1)' shapes.BoundingBox(:,2)'];
        end
    else
        fprintf(err_msg);
        return;
    end
end
if size(centroids_rect) ~= 4,       fprintf(err_msg); return; end
if ~isnumeric(centroids_rect),      fprintf(err_msg); return; end
if ~exist('resolution_km','var') || isempty(resolution_km)
    fprintf('No resolution specified, resorting to default 1 km \n');
    resolution_km = 1;
end
resolution_ang = resolution_km / (111.12);

min_lon = centroids_rect(1);
max_lon = centroids_rect(2);
min_lat = centroids_rect(3);
max_lat = centroids_rect(4);

n_lon = round((max_lon - min_lon)/(resolution_ang)) + 1;
n_lat = round((max_lat - min_lat)/(resolution_ang)) + 1;
n_centroids = n_lon * n_lat;

fprintf(sprintf('generating centroids at %3.2f km resolution \n', resolution_km));
for i = 0 : n_lon - 1
    ndx = i * n_lat;
    centroids.lat(1,ndx + 1 : ndx + n_lat)= (1:n_lat) .* resolution_ang + min_lat;
    centroids.lon(1,ndx + 1 : ndx + n_lat)= (n_lon - i) .* resolution_ang + min_lon;
end
centroids.centroid_ID = [1:n_centroids];

centroids.onLand = true(size(centroids.centroid_ID));

if exist(climada_global.map_border_file,'file')
    load(climada_global.map_border_file)
    in = inpolygon(centroids.lon,centroids.lat,shapes.X,shapes.Y);
    centroids.onLand = false(size(centroids.centroid_ID));
    centroids.onLand(in) = 1;

    % centroids.onLand = true(size(centroids.centroid_ID));
    
    if isfield(shapes,'NAME')
        for i = 1 : n_centroids
            centroids.countryname{i} = shapes.NAME;
        end
        centroids.admin0_name = shapes.NAME;
    end
    if isfield(shapes,'ADM0_A3')
        centroids.admin0_ISO3 = shapes.ADM0_A3;
    end
else
    fprintf('WARNING: no border info found, centroids.onLand set to 1 for all centroids \n')
end

