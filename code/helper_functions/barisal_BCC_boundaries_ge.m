function barisal_BCC_boundaries_ge(google_earth_save)
% climada
% MODULE:
%   barisal_demo
% NAME:
%   barisal_BCC_boundaries_ge
% PURPOSE:
%   visualisation of BCC boundaries in google earth
% CALLING SEQUENCE:
%   barisal_BCC_boundaries_ge(google_earth_save)
% EXAMPLE:
%   barisal_BCC_boundaries_ge
% INPUTS:
%   google_earth_save: the filename of the resulting .kmz google earth file
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150315, initial
%-

% res=[]; % init output
close all % not really necessary, but speeds things up

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('google_earth_save','var'), google_earth_save = []; end

% PARAMETERS
google_earth_save = 'BCC_boundary';


%% read BCC admin boundary shapefile (WGS84 projection) and save as mat-file
% shp_file     = 'M:\BGCC\CHR\RK\RS\A_Sustainable_Development\Projects\ECA\BarisalBangladesh\Barisal_GIS\WGS1984\AdminBoundaryLineBCC_Project.shp';
% BCC_boundary = climada_shaperead(shp_file,0,1,0,1); % 0 for no-save
% shp_mat_file = [climada_global.data_dir filesep 'results' filesep 'BCC_boundary_shp_in_GIS.mat'];
% save(shp_mat_file, 'BCC_boundary')


%% read BCC admin boundary shapefile (original GCS Everest 1830 coordinates, UTM)
shp_file     = 'M:\BGCC\CHR\RK\RS\A_Sustainable_Development\Projects\ECA\BarisalBangladesh\Barisal_GIS\GIS\Admin Boundary Line BCC.shp';
BCC_boundary = climada_shaperead(shp_file,0,1,0,1); % 0 for no-save


% copy X to X_BTM and Y to Y_BTM (to keep original data)
% and convert to lat lon (incl. manual shift)
for shape_i = 1:length(BCC_boundary)
    BCC_boundary(shape_i).X_BTM = BCC_boundary(shape_i).X;
    BCC_boundary(shape_i).Y_BTM = BCC_boundary(shape_i).Y;
    [BCC_boundary(shape_i).Y, BCC_boundary(shape_i).X] = utm2ll_shift(BCC_boundary(shape_i).X_BTM, BCC_boundary(shape_i).Y_BTM);
end
shp_mat_file = [climada_global.data_dir filesep 'results' filesep 'BCC_boundary_shp.mat'];
save(shp_mat_file, 'BCC_boundary')



%% prompt for google_earth_save file
if isempty(google_earth_save) % local GUI
    google_earth_save = [climada_global.data_dir filesep 'results' filesep 'Select name to save google earth visualiation .kmz'];
    [filename, pathname] = uiputfile(google_earth_save, 'Save google earth as:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        google_earth_save = fullfile(pathname,filename);
    end
end

[pathstr, name, ext] = fileparts(google_earth_save);
if ~strcmp(ext,'kmz')
    ext = '.kmz';
end
if strcmp(pathstr,'')
    pathstr = [climada_global.data_dir filesep 'results'];
end
google_earth_save = [fullfile(pathstr,name) ext];

% save kmz file
fprintf('saving BCC boundary as\n %s\n',google_earth_save); 
k = kml(google_earth_save);


%% create google earth output
kk = k.newFolder('BCC boundaries');
for shape_i=1:length(BCC_boundary)
    X = BCC_boundary(shape_i).X(~isnan(BCC_boundary(shape_i).X));
    Y = BCC_boundary(shape_i).Y(~isnan(BCC_boundary(shape_i).Y));
    kk.plot(X,Y); %'lineColor','#00FFFFFF'
end % shape_i


%% open visualition in google earth
k.run  



    
    
    

