function climada_shp_explorer(shapes_folder,single_mode,fast_check)
% climada
% MODULE
%   _shapes
% NAME:
%   climada_shp_explorer
% PURPOSE:
%   explore content of shape files and folders full of shape files
%
%   WARNING: might take a LONG time, hence see fast_check
%
% CALLING SEQUENCE:
%   climada_shp_explorer(shapes_folder,single_mode,fast_check)
% EXAMPLE:
%   climada_shp_explorer('',0,10) % plot every 10th shape (often a good start)
% INPUTS:
%   shapes_folder: a folder with shape files (or a single shape file, if
%       extension .shp, if single_mode=1)
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   single_mode: =0, select all files in the folder
%       =1 only the actually selected file
%   fast_check: >1: only plot every fast check shape (to get a feel)
%       start e.g. with fast_check=100 and then fast_check=10
%       =1: plot all shapes (default)
% OUTPUTS:
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20141223
%-

%global climada_global % access to global variables
if ~climada_init_vars,return;end % init/import global variables

%%if climada_global.verbose_mode,fprintf('*** %s ***\n',mfilename);end % show routine name on stdout

% poor man's version to check arguments
if ~exist('shapes_folder','var'),shapes_folder='';end
if ~exist('single_mode','var'),single_mode=0;end
if ~exist('fast_check','var'),fast_check=1;end

% locate the module's data
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% PARAMETERS
%
% set list of colors for lines
color_list='rgbcmyk';
color_list_long={'Black','Cyan','Magenta','Blue','Green','Red','Yellow'};

% prompt for shapes_folder if not given
if isempty(shapes_folder) % local GUI
    shapes_folder=[module_data_dir filesep '*.shp'];
    [filename, pathname] = uigetfile(shapes_folder, 'Select first shape file:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        shapes_folder=fullfile(pathname,filename);
    end
end

if single_mode
    shapes_files{1}=shapes_folder;
else
    % create list of shape files
    shapes_files={};
    fP=fileparts(shapes_folder);
    D=dir([fP filesep '*.shp']);
    for D_i=1:length(D)
        if ~D(D_i).isdir
            shapes_files{end+1}=[fP filesep D(D_i).name];
        end
    end % D_i
end % single_mode

fprintf('processing %i files:\n',length(shapes_files))

% plotting frenzy starts
% ----------------------

color_i=1; % init
for file_i=1:length(shapes_files)
    [~,fN]=fileparts(shapes_files{file_i});
    fprintf('adding %s ',fN)
    shapes=climada_shaperead(shapes_files{file_i},0,0,0,1); % 0 for no-save
    fprintf('(#%i, %s): ',length(shapes),color_list(color_i))
    field_names = fieldnames(shapes(1));
    for field_i=1:length(field_names)
        fprintf('%s ',field_names{field_i})
    end
    fprintf('\n')
    if strcmp(fN,'natural')
        color_i=color_i+1;
        if color_i>length(color_list),color_i=1;end
        continue;
    end
    for shape_i=1:fast_check:length(shapes)
        if strcmp(shapes(shape_i).Geometry,'Polygon')
            pos=find(~isnan(shapes(shape_i).X)); % remove NaN to fill
            fill(shapes(shape_i).X(pos),shapes(shape_i).Y(pos),color_list(color_i)),hold on
        else % most others either 'Point' or 'Line'
            if strcmp(fN,'buildings')
                plot3(shapes(shape_i).X,shapes(shape_i).Y,zeros(size(shapes(shape_i).X))+10, ['s' color_list(color_i)]),hold on  
            else
                plot3(shapes(shape_i).X,shapes(shape_i).Y,zeros(size(shapes(shape_i).X))+10, ['-' color_list(color_i)]),hold on
            end
        end
%         if strcmp(shapes(shape_i).Geometry,'Point') % label points, but often TOO MESSY
%             text(shapes(shape_i).X,shapes(shape_i).Y,shapes(shape_i).name);
%         end
    end % shape_i
    color_i=color_i+1;
    if color_i>length(color_list),color_i=1;end
end % file_i
axis equal

% climada_plot_world_borders(2,'','',1)
set(gcf,'Color',[1 1 1])
drawnow

end
