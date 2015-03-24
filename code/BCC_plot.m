function BCC_plot
% climada
% MODULE
%   barisal_demo
% NAME:
%   Barisal_plot
% PURPOSE:
%   plot Barisal map, including all relevant information (BCC GIS and open
%   data). Hardwired for Barisal.
% CALLING SEQUENCE:
%   Barisal_plot
% EXAMPLE:
%   Barisal_plot
% INPUTS:
%   
% OPTIONAL INPUT PARAMETERS:
%   
% OUTPUTS:
% MODIFICATION HISTORY:
% Lea Mueller, 20150305, muellele@gmail.com
%-

global climada_global % access to global variables
if ~climada_init_vars,return;end % init/import global variables


%parameters
check_printplot = 0;
axislim = [90.25 90.45 22.6 22.8]; %barisal close up BCC
% axislim = [90.20 90.55 22.55 22.85]; %barisal bigger BCC
% axislim = [88.5 93 21 24]; %coastal bangladesh


% locate the module's data
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
% for testing
% module_data_dir  = '\\CHRB1065.CORP.GWPNET.COM\homes\X\S3BXXW\Documents\lea\climada_git\climada_modules\barisal_demo\data';
GIS_dir          = 'M:\BGCC\CHR\RK\RS\A_Sustainable_Development\Projects\ECA\BarisalBangladesh\Barisal_GIS\WGS1984';
GIS_open_dat_dir = [module_data_dir filesep 'entities' filesep 'BGD_adm'];


% read shape files
% BCC boundaries
BCC_file = 'CityCorporationAreaPolyBCC_P.shp'; 
BCC      = climada_shaperead([GIS_dir filesep BCC_file],0,1,0,1); % 0 for no-save

% admin 4 
shp_file_ = 'BGD_adm'; i = 4;
shp_file = sprintf('%s%d.shp',shp_file_,i);
shapes   = climada_shaperead([GIS_open_dat_dir filesep shp_file],0,1,0,1); % 0 for no-save
    

% load 
BCC_entity_file = [module_data_dir filesep 'entities' filesep 'Barisal_BCC_1km_100.mat'];
if exist(BCC_entity_file,'file')
    load(BCC_entity_file)
end



%% prepare figure
close all; h=[];
fig = climada_figuresize(0.5,0.7);
h(1) = plot(BCC(1).X,BCC(1).Y,'color',[240 128 128]/255);
hold on
% admin4
for shape_i=1:length(shapes)
    h(2)= plot(shapes(shape_i).X,shapes(shape_i).Y,'color',[191 191 191]/255);%grey
end
climada_plot_world_borders(0.5)
axis equal
axis(axislim)
title('Barisal')
legend(h,'BCC boundaries','admin4')
climada_plot_entity_assets(entity,'','','','',1);


printname = '';
if check_printplot
    if isempty(printname) %local GUI
        printname_         = [climada_global.data_dir filesep 'results' filesep '*.pdf'];
        printname_default  = [climada_global.data_dir filesep 'results' filesep 'type_name.pdf'];
        [filename, pathname] = uiputfile(printname_,  'Save asset map as figure:',printname_default);
        foldername = [pathname filename];
        if pathname <= 0; return;end
    else
        foldername = [climada_global.data_dir filesep 'results' filesep 'Entity.pdf'];
    end
    print(fig,'-dpdf',foldername)
    cprintf([255 127 36 ]/255,'\t\t saved 1 FIGURE in folder ..%s \n', foldername);
end
    
    

return



%% for all admin layers

% colors_ = jet(4);
% % admin1,2,3,4
% shp_file_ = 'BGD_adm';
% for i = 4:-1:1
%     shp_file = sprintf('%s%d.shp',shp_file_,i);
%     shapes   = climada_shaperead([GIS_open_dat_dir filesep shp_file],0,1,0,1); % 0 for no-save
%     for shape_i=1:length(shapes)
%         if strcmp(shapes(shape_i).Geometry,'Polygon')
%             h(i+1) = plot(shapes(shape_i).X,shapes(shape_i).Y,'color',colors_(i,:));
%             
%             %pos=find(~isnan(shapes(shape_i).X)); % remove NaN to fill
%             %plot(shapes(shape_i).X(pos),shapes(shape_i).Y(pos));
%             %fill(shapes(shape_i).X(pos),shapes(shape_i).Y(pos),color_list(color_i)),hold on
%         else % most others either 'Point' or 'Line'
%             %if strcmp(fN,'buildings')
%             %    plot3(shapes(shape_i).X,shapes(shape_i).Y,zeros(size(shapes(shape_i).X))+10, ['s' color_list(color_i)]),hold on  
%             %else
%             %    plot3(shapes(shape_i).X,shapes(shape_i).Y,zeros(size(shapes(shape_i).X))+10, ['-' color_list(color_i)]),hold on
%             %end
%         end
%     end % shape_i
% end
% 
% legend(h,'BCC boundaries','admin4', 'admin3', 'admin2', 'admin1')


