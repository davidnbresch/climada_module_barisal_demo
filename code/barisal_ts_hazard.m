
%% hazard tc wind (historical): just loading not calculating
hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC'];
load(hazard_set_file)

% hazard tropical cyclone surge (ts)
hazard_set_file_ts = strrep(hazard_set_file,'_TC', '_TS.mat');


%% the shape file with higher resolution for Bangladesh
module_data_dir = [climada_global.modules_dir filesep 'barisal_demo' filesep 'data']; 
BGD_country_shapefile=[module_data_dir filesep 'entities' filesep 'BGD_adm' filesep 'BGD_adm0.shp'];
[fP,fN] = fileparts(BGD_country_shapefile);
BGD_country_shapefile_mat=[fP filesep fN '.mat'];
if ~exist(BGD_country_shapefile_mat,'file')
    % first time, read the shape file and store as .mat
    climada_shaperead(BGD_country_shapefile);
end

% replace the map border shape file (be careful, they should be reset)
climada_global.map_border_file = BGD_country_shapefile_mat;

load(BGD_country_shapefile_mat)

%% DEM
adm_lvl = 4;
DEM_save_file = [module_data_dir filesep 'system' filesep strcat('Barisal_',num2str(adm_lvl),'_DEM.mat')];
load(DEM_save_file) 

% DEM figure
% climada_DEM_plot()
figure
hold on
for shape_i=1:length(shapes)
    h(2)= plot(shapes(shape_i).X,shapes(shape_i).Y,'color',[191 191 191]/255);%grey
end
centroids.lon = DEM.lon;
centroids.lat = DEM.lat;
[X, Y, gridded_VALUE] = climada_gridded_VALUE(DEM.elevation_m,centroids);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
colorbar

    

%% --centroids-----------
% load centroids
centroids_file  = [climada_global.data_dir filesep 'system' filesep 'Barisal_BCC_centroids'];
load(centroids_file)
% figure
% plot(centroids.lon, centroids.lat,'+')

%% assign elevation to centroids
% srtm_data_dir = [module_data_dir filesep 'system' filesep 'srtm_55_08'];
% check_plots = 1;
% [DEM, centroids] = climada_read_srtm_DEM(srtm_data_dir,centroids, DEM_save_file, 5, check_plots);
% [DEM, centroids] = climada_read_srtm_DEM(srtm_dir, centroidsORcountry, DEM_save_file, smooth, check_plot)

% hazard_tc.elevation_m = interp2(DEM.lon,DEM.lat,DEM.elevation_m,hazard_tc.lon(1),hazard_tc.lat(1),'nearest');

hazard.elevation_m = zeros(size(hazard.lon));
cos_tc_track_lat   = cos(DEM.lat/180*pi);
for i= 1:numel(hazard.lon)
    dd=((DEM.lon-hazard.lon(i)).*cos_tc_track_lat).^2+(DEM.lat-hazard.lat(i)).^2; % in km^2
    [~,pos] = min(dd);
    node_i  = pos(1); % take first if more than one
    D = sqrt(dd(node_i))*111.12; % now in km
    hazard.elevation_m(i) = DEM.elevation_m(node_i);
end
    
% check hazard elevation data with figure
centroids_short = centroids;
centroids_short.lon = centroids_short.lon(1:64);
centroids_short.lat = centroids_short.lat(1:64);
figure
[X, Y, gridded_VALUE] = climada_gridded_VALUE(hazard.elevation_m(1:64),centroids_short);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
colorbar
hold on
for shape_i=1:length(shapes)
    h(2)= plot(shapes(shape_i).X,shapes(shape_i).Y,'color',[191 191 191]/255);%grey
end
plot(centroids_short.lon, centroids_short.lat, '+')


%% calculate ts hazard
hazard.onLand = ones(size(hazard.lon));
hazard_tc = hazard;
hazard_ts = tc_surge_hazard_create(hazard_tc,hazard_set_file_ts,DEM);


% plot sidr
event_i = 173;
label = ''; 
caxis_range = [0 8];
plot_centroids = 1;
% figure
% res=climada_hazard_plot(hazard_ts,event_i,label,caxis_range,plot_centroids)
figure
[X, Y, gridded_VALUE] = climada_gridded_VALUE(full(hazard_ts.intensity(event_i,1:64)),centroids_short);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
for shape_i=1:length(shapes)
    h(2)= plot(shapes(shape_i).X,shapes(shape_i).Y,'color',[191 191 191]/255);%grey
end
colorbar
plot(centroids_short.lon, centroids_short.lat, '+')
plot(centroids.lon(65:100), centroids.lat(65:100), 'r+')
plot(centroids.lon(101:200), centroids.lat(101:200), 'r+')


%% define centroids in the river --> centroids_ts
climada_define_polygon
centroids_ts = centroids;
centroids_ts.lon = [centroids.lon(1:64) ans(:,1)'];
centroids_ts.lat = [centroids.lat(1:64) ans(:,2)'];
centroids_ts.comment = 'handmade with river, 24-Mar-2015';
centroids_ts.centroid_ID = 1:numel(centroids_ts.lon);
centroids_file_ts  = [centroids_file '_ts'];
save(centroids_file_ts,'centroids_ts');

% calculate tc hazard
hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_TC_centroids_TS'];
hazard = climada_tc_hazard_set(tc_track, hazard_set_file, centroids_ts);

% add elevation information to hazard and save tc hazard
hazard.elevation_m = zeros(size(hazard.lon));
cos_tc_track_lat   = cos(DEM.lat/180*pi);
for i= 1:numel(hazard.lon)
    dd=((DEM.lon-hazard.lon(i)).*cos_tc_track_lat).^2+(DEM.lat-hazard.lat(i)).^2; % in km^2
    [~,pos] = min(dd);
    node_i  = pos(1); % take first if more than one
    D = sqrt(dd(node_i))*111.12; % now in km
    hazard.elevation_m(i) = DEM.elevation_m(node_i);
end
save(hazard_set_file,'hazard')

% plot elevation
figure
[X, Y, gridded_VALUE] = climada_gridded_VALUE(hazard.elevation_m,centroids_ts);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
colorbar
hold on
for shape_i=1:length(shapes)
    h(2)= plot(shapes(shape_i).X,shapes(shape_i).Y,'color',[191 191 191]/255);%grey
end
plot(centroids_ts.lon, centroids_ts.lat, '+')
c_i = 12;
plot(centroids_ts.lon(64+c_i), centroids_ts.lat(64+c_i), 'rd','markersize',10)

% calculate ts hazard
hazard.onLand = ones(size(hazard.lon));
hazard_tc = hazard;
hazard_ts = tc_surge_hazard_create(hazard_tc,hazard_set_file_ts,DEM);

% plot sidr
event_i = 173;
label = ''; 
caxis_range = [0 8];
plot_centroids = 1;
% figure
% res=climada_hazard_plot(hazard_ts,event_i,label,caxis_range,plot_centroids)
figure
[X, Y, gridded_VALUE] = climada_gridded_VALUE(full(hazard_ts.intensity(event_i,:)),centroids_ts);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
for shape_i=1:length(shapes)
    h(2)= plot(shapes(shape_i).X,shapes(shape_i).Y,'color',[191 191 191]/255);%grey
end
colorbar
plot(centroids_ts.lon, centroids_ts.lat, '+')


%% plot 1998 cyclone, track 127
event_i = 127;
event_i = 173;
label = ''; 
caxis_range = [0 8];
plot_centroids = 1;
% figure
% res=climada_hazard_plot(hazard_ts,event_i,label,caxis_range,plot_centroids)
figure
[X, Y, gridded_VALUE] = climada_gridded_VALUE(full(hazard_ts.intensity(event_i,:)),centroids_ts);
contourf(X, Y, gridded_VALUE,200,'edgecolor','none')
hold on
for shape_i=1:length(shapes)
    h(2)= plot(shapes(shape_i).X,shapes(shape_i).Y,'color',[191 191 191]/255);%grey
end
axis([90.30 90.4 22.64 22.72])
axis equal
colorbar
plot(centroids_ts.lon, centroids_ts.lat, '+')

% create IFC for surge heights at centroid 12+64
c_i = 12;
hazard_ts.elevation_m(64+c_i)
IFC = climada_hazard2IFC(hazard_ts,64+c_i);
figure
climada_IFC_plot(IFC)

for c_i = 1:20 %-63:0%20:30%
    event_i = 127;
    a = full(hazard_ts.intensity(event_i,64+c_i));
    fprintf('%d: 1998 surge %2.2f\n', c_i, a);
    event_i = 173;
    a = full(hazard_ts.intensity(event_i,64+c_i));
    fprintf('%d: 2007 surge %2.2f\n', c_i, a);
end



%% print a figure
% foldername  = [filesep 'results' filesep 'Sidr_ts_footprint_ts_additional_centroids.pdf'];
% print(gcf,'-dpdf',[climada_global.data_dir foldername])



