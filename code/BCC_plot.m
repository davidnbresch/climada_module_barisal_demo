function BCC_plot
% climada
% MODULE
%   barisal_demo
% NAME:
%   BCC_plot
% PURPOSE:
%   plot Barisal map, including all relevant information (BCC GIS, BCC wards 
%   polygon and open data). Hardwired for Barisal.
% CALLING SEQUENCE:
%   BCC_plot
% EXAMPLE:
%   BCC_plot
% INPUTS:
%   none, all hardwired for Barisal
% OUTPUTS:
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150305, 
% Lea Mueller, muellele@gmail.com, 20150506, add new BCC ward polygons
%-

global climada_global % access to global variables
if ~climada_init_vars,return;end % init/import global variables


%parameters
% check_printplot = 0;
check_printplot = 1;
axislim = [90.28 90.41 22.64 22.775]; %barisal close up BCC
% axislim = [90.285 90.3957 22.64 22.752]; %barisal close up BCC 
% axislim = [90.25 90.45 22.6 22.8]; %barisal close up BCC
% axislim = [90.20 90.55 22.55 22.85]; %barisal bigger BCC
% axislim = [88.5 93 21 24]; %coastal bangladesh
axis(axislim)

% locate the module's data
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
% for testing
% module_data_dir  = '\\CHRB1065.CORP.GWPNET.COM\homes\X\S3BXXW\Documents\lea\climada_git\climada_modules\barisal_demo\data';
% GIS_dir          = 'M:\BGCC\CHR\RK\RS\A_Sustainable_Development\Projects\ECA\BarisalBangladesh\Barisal_GIS\WGS1984';
GIS_dir          = 'M:\BGCC\CHR\RK\RS\A_Sustainable_Development\Projects\ECA\BarisalBangladesh\Barisal_GIS\GIS';
GIS_open_dat_dir = [module_data_dir filesep 'entities' filesep 'BGD_adm'];


%% read shape files
%% BCC boundaries
% BCC_file = 'CityCorporationAreaPolyBCC_P.shp'; 
BCC_file = 'City Corporation Area Poly BCC.shp'; 
BCC      = climada_shaperead([GIS_dir filesep BCC_file],0,1,0,1); % 0 for no-save
%transform local coordinates (GCS  Everest 1830) to WGS1984
[BCC(1).lon,BCC(1).lat] = utm2ll_shift(BCC(1).X, BCC(1).Y);
%save BCC
BCC_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_border.mat'];
save(BCC_savename,'BCC')


%% BCC Wards
% ward_file = 'Ward Boundary Poly BCC.shp';
% ward_file = 'Ward Boundary Line BCC.shp';
% ward_file = 'Wards_popdens.shp';
ward_file = 'Ward_popdens.shp';
BCC_wards = climada_shaperead([GIS_dir filesep ward_file],0,1,0,1); % 0 for no-save
%transform local coordinates (GCS  Everest 1830) to WGS1984
for w_i = 1:length(BCC_wards)
    [BCC_wards(w_i).lon,BCC_wards(w_i).lat] = utm2ll_shift(BCC_wards(w_i).X, BCC_wards(w_i).Y);
end
for w_i = 1:length(BCC_wards)
    BCC_wards(w_i).lon(isnan(BCC_wards(w_i).lon)) = [];
    BCC_wards(w_i).lat(isnan(BCC_wards(w_i).lat)) = [];
end
%save BCC_wards
BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards.mat'];
save(BCC_wards_savename,'BCC_wards')
BCC_wards_temp = BCC_wards;

% % add ward numbers from previous file - not possible, the ordering is different!
% BCC_wards_savename2 = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards_number_added.mat'];
% load(BCC_wards_savename2)
% for w_i = 1:length(BCC_wards_temp)
%     BCC_wards_temp(w_i).Ward_no = BCC_wards(w_i).Ward_no;
% end
% BCC_wards = BCC_wards_temp;
% save(BCC_wards_savename,'BCC_wards')



%% BCC Wards 2
% ward_file   = ['Barisal_ATLAS - Excel and GIS' filesep 'Barisal' filesep 'All' filesep 'Ward_boundary.shp'];
% folder_name = ['Vulnerability_Analysis_Barisal' filesep '2_data' filesep 'reports_and_GIS_files']; 
% BCC_wards   = climada_shaperead([strrep(climada_global.data_dir,[filesep 'climada_data'],'') filesep folder_name filesep ward_file],0,1,0,1); % 0 for no-save
% %transform local coordinates (GCS  Everest 1830) to WGS1984
% for w_i = 1:length(BCC_wards)
%     [BCC_wards(w_i).lon,BCC_wards(w_i).lat] = utm2ll_shift(BCC_wards(w_i).X, BCC_wards(w_i).Y);
% end
% for w_i = 1:length(BCC_wards)
%     BCC_wards(w_i).lon(isnan(BCC_wards(w_i).lon)) = [];
%     BCC_wards(w_i).lat(isnan(BCC_wards(w_i).lat)) = [];
% end


%% admin 4, open GIS data
shp_file_= 'BGD_adm'; i = 4;
shp_file = sprintf('%s%d.shp',shp_file_,i);
shapes   = climada_shaperead([GIS_open_dat_dir filesep shp_file],0,1,0,1); % 0 for no-save
    

%% load generic 10km entity
BCC_entity_file = [module_data_dir filesep 'entities' filesep 'Barisal_BCC_1km_100.mat'];
if exist(BCC_entity_file,'file')
    load(BCC_entity_file)
end



%% prepare figure
close all; h=[];
fig = climada_figuresize(0.5,0.7);
% h(1) = plot(BCC(1).X,BCC(1).Y,'color',[240 128 128]/255);
h(1) = plot(BCC(1).lon,BCC(1).lat,'color',[240 128 128]/255);
hold on

% plot admin4
for shape_i=1:length(shapes)
    h(2)= plot(shapes(shape_i).X,shapes(shape_i).Y,'color',[191 191 191]/255);%grey
end

% plot BCC wards (fill polygons)
color_ = jet(length(BCC_wards));
for w_i=1:length(BCC_wards)
    %h(3)= plot(BCC_wards(w_i).lon,BCC_wards(w_i).lat,'color',[244 164 96 ]/255,'linewidth',2);%sandybrown
    h(3)= plot(BCC_wards(w_i).lon,BCC_wards(w_i).lat,'color','k','linewidth',1);%sandybrown
    %fill(BCC_wards(w_i).lon,BCC_wards(w_i).lat,color_(w_i,:));%sandybrown
end
% add labels (ward names)
for w_i=1:length(BCC_wards)
    text(mean(BCC_wards(w_i).lon), mean(BCC_wards(w_i).lat), BCC_wards(w_i).UNION_NAME)
end

climada_plot_world_borders(0.5)
axis equal
axis(axislim)
title('Barisal')
legend(h,'BCC boundaries','admin4','BCC wards')
% climada_plot_entity_assets(entity,'','','','',1);

printname = '';
if check_printplot
    if isempty(printname) %local GUI
        printname_         = [climada_global.data_dir filesep 'results' filesep '*.pdf'];
        printname_default  = [climada_global.data_dir filesep 'results' filesep 'type_name.pdf'];
        [filename, pathname] = uiputfile(printname_,  'Save BCC figure as:',printname_default);
        foldername = [pathname filename];
        if pathname <= 0; return;end
    else
        foldername = [climada_global.data_dir filesep 'results' filesep 'BCC_ward_plot.pdf'];
    end
    print(fig,'-dpdf',foldername)
    cprintf([255 127 36 ]/255,'\t\t saved 1 FIGURE in folder ..%s \n', foldername);
end
    

% % ward
% for w_i=1:length(BCC_wards2)
%     h(3)= plot(BCC_wards2(w_i).lon,BCC_wards2(w_i).lat,'color',[244 164 96 ]/255,'linewidth',2);%sandybrown
%     %h(3)= fill(BCC_wards(w_i).lon,BCC_wards(w_i).lat,[244 164 96 ]/255);%sandybrown
% end
% figure
% w_i= 1;
% h(3)= plot(BCC_wards(w_i).lon,BCC_wards(w_i).lat,'color','b','linewidth',2);%sandybrown
% hold on
% fill(BCC_wards(w_i).X,BCC_wards(w_i).Y,[244 164 96 ]/255);%sandybrown

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


