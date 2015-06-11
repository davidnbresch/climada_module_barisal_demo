function fig=climada_EDS_plot_3d(hazard,EDS,event_i,check_labels,varargin)
% climada visualise event damage
% NAME:
%   climada_plot_EDS_3d
% PURPOSE:
%   Plot the damage to assets from event damage set (EDS - calculated using
%   climada_EDS_calc) as a 3d bar chart, overlaying the country borders.
%   Further, a surface (for TS) or a contour (for other hazards) is plotted
%   to show the spatial extent.
% CALLING SEQUENCE:
%   fig = climada_plot_EDS_3d(hazard,EDS,event_i)
% EXAMPLE:
%   climada_plot_EDS_3d(hazard,EDS,event_i)
% INPUTS:
%   hazard: hazard structure
%   EDS: event damage set structure, see climada_EDS_calc
% OPTIONAL INPUT PARAMETERS:
%   event_i: specify event number of interest (e.g. use climada_hazard_plot
%       to figure ID of largest, 2nd largest... event)
%       If empty, the largest event is displayed
% OUTPUTS:
%   fig             : figure handle
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com 20141204
% David N. Bresch, david.bresch@gmail.com, 20141215, cleanup
% Gilles Stassen, 20150114 fixed conversion issues + cleanup
% Gilles Stassen, 20150517 added varargin for plotting extra features, fixed issue with axis labeling
%-
fig = []; % init output

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('hazard'  ,       'var'), return;         end
if ~exist('EDS'     ,       'var'), return;         end
if ~exist('event_i' ,       'var'), event_i = -1;   end % Default: pick largest event
if ~exist('check_labels',   'var'), check_labels=1; end
% Prepare data
if event_i<0
    % search for i-th largest event
    event_sum=sum(hazard.intensity,2);
    [~,sorted_i]=sort(event_sum);
    event_ii=sorted_i(length(sorted_i)+event_i+1);
    hazard_intensity=full(hazard.intensity(event_ii,:)); % extract one event
    if event_i<-1
        title_str=sprintf('%s %i-largest event (%i)',hazard.peril_ID,-event_i,event_ii);
    else
        title_str=sprintf('%s largest event (%i)',hazard.peril_ID,event_ii);
    end
    % plot some further info to sdout:
    if (isfield(hazard,'name') && isfield(hazard,'yyyy')) && (isfield(hazard,'mm') && isfield(hazard,'dd'))
        fprintf('%s, %4.4i%2.2i%2.2i, event %i\n',hazard.name{event_ii},hazard.yyyy(event_ii),hazard.mm(event_ii),hazard.dd(event_ii),event_ii);
    end
    event_i = event_ii;
elseif event_i==0
    hazard_intensity=full(max(hazard.intensity)); % max intensity at each point
    title_str=sprintf('%s max damage at each centroid',hazard.peril_ID);
else
    hazard_intensity=full(hazard.intensity(event_i,:)); % extract one event
    title_str=sprintf('%s event %i',hazard.peril_ID,event_i);
    % plot some further info to sdout:
    if (isfield(hazard,'name') && isfield(hazard,'yyyy')) && (isfield(hazard,'mm') && isfield(hazard,'dd'))
        fprintf('%s, %4.4i%2.2i%2.2i, event %i\n',hazard.name{event_i},hazard.yyyy(event_i),hazard.mm(event_i),hazard.dd(event_i),event_i);
    end
end

% Matlab's bar3 function takes a matrix and plots it without being able to
% define x and y coords. The following is a round about method of scaling
% all order parts of the plot to the axes defined by bar3
% -------------------

% First, generate matrix from singleton asset arrays. Define the x and y
% axes as tmp_lon and tmp_lat respectively
tmp_lon = unique(EDS.assets.lon);
tmp_lat = unique(EDS.assets.lat);

% tmp_lon = [min(EDS.assets.lon):abs(min(diff(EDS.assets.lon))):max(EDS.assets.lon)];
% tmp_lat = [min(EDS.assets.lat):abs(min(diff(EDS.assets.lat))):max(EDS.assets.lat)];

% Conversion factors
dx_dlon = numel(tmp_lon)/(max(tmp_lon) - min(tmp_lon));
dy_dlat = numel(tmp_lat)/(max(tmp_lat) - min(tmp_lat));

lon_c =  min(tmp_lon);  lat_c = min(tmp_lat);

% Construct matrix
loss_plot = zeros(numel(tmp_lat),numel(tmp_lon));
if event_i == 0
    tmp_d_a_c = EDS.ED_at_centroid;
else
    tmp_d_a_c = full(EDS.damage_at_centroid(:,event_i));
end
x = zeros(size(EDS.assets.Value));
y = zeros(size(EDS.assets.Value));

for ndx = 1: numel(EDS.assets.Value)
    x = find(tmp_lon == EDS.assets.lon(ndx));
    y = find(tmp_lat == EDS.assets.lat(ndx));
    loss_plot(y,x) = tmp_d_a_c(ndx);
end

% Plot 3d bar chart
fig = figure('Name',title_str,'Color',[1 1 1]);
set(fig,'name', title_str);
hold on
if event_i == 0
    title(['Annual expected damage ' num2str(EDS.reference_year) ': USD '  num2str(round(EDS.ED/1e3)) 'k'],'fontsize',12)
else
    title(['Total damage for event ' num2str(event_i) ': USD '  num2str(round(sum(EDS.damage_at_centroid(:,event_i))/1e3)) 'k'],'fontsize',12)
end
%title(title_str)
h = bar3(loss_plot);
for i = 1:size(loss_plot,2)
    z = get(h(i),'z');
    k = 1;
    for j = 0:6:(6*size(loss_plot,1)-6)
        if loss_plot(k,i)<=0
            z(j+1:j+6,:)= NaN; % set zeroes to NaN to only show non-zero bars
        end
        k = k+1;
    end
    set(h(i),'z',z);
end

% Set color data equal to height data to color bars based on value
for k = 1:length(h)
    zdata = get(h(k),'zdata');
    set(h(k),'cdata',zdata);
    set(h(k),'FaceColor','interp');
end

% Color bar form light (ish) blue to red
R = linspace(0,1,100);
G = linspace(0.5,0,100);
B = linspace(1,0,100);

map = [R' G' B'];

colormap(climada_colormap('schematic'));

hold on;

% Plot borders
load(climada_global.map_border_file)

% store also in one contiguous list (wrap-around for backward
% compatibility)
whole_world_borders.lon = [NaN];
whole_world_borders.lat = [NaN];
for i=1:length(shapes)
    whole_world_borders.lon = [whole_world_borders.lon; NaN; shapes(i).X']; % NaN at end already there
    whole_world_borders.lat = [whole_world_borders.lat; NaN; shapes(i).Y']; % NaN at end already there
end

% World borders to area of interest
tmp_lat_wb = (whole_world_borders.lat >= min(tmp_lat) &...
    whole_world_borders.lat <= max(tmp_lat)) | isnan(whole_world_borders.lat);
tmp_lon_wb = (whole_world_borders.lon >= min(tmp_lon) &...
    whole_world_borders.lon <= max(tmp_lon)) | isnan(whole_world_borders.lon);

whole_world_borders.lat = whole_world_borders.lat(tmp_lat_wb & tmp_lon_wb);
whole_world_borders.lon = whole_world_borders.lon(tmp_lat_wb & tmp_lon_wb);

% Scale borders to match bar3
x_b = dx_dlon.*(whole_world_borders.lon - lon_c);
y_b = dy_dlat.*(whole_world_borders.lat - lat_c);

plot(x_b,y_b,'color',[81 81 81]/256);

if ~isempty(varargin),      shape_file_i = 0;   s = [];      end
for arg_i = 1:4:length(varargin)
    err_msg = sprintf('ERROR: invalid argument\n');
    if ~ischar(varargin{arg_i}),        cprintf([1 0 0],err_msg);  return;     end
    if length(varargin{arg_i+1}) ~= 1,  cprintf([1 0 0],err_msg);  return;     end
    if length(varargin{arg_i+2}) ~= 3,  cprintf([1 0 0],err_msg);  return;     end
    if ~ischar(varargin{arg_i+3}),      cprintf([1 0 0],err_msg);  return;     end

    [fP,fN,fE] = fileparts(varargin{arg_i});
    if strcmp(fE,'.shp')
        shapes = shaperead([fP filesep fN fE]);
    elseif strcmp(fE,'.mat')
        if exist('shapes','var'),   clear shapes;   end
        load([fP filesep fN fE])
        if ~exist('shapes','var')
            cprintf([1 0 0], 'ERROR: data in .mat file must have variable name ''shapes''\n');
            return
        end
    else
        cprintf([1 0 0],err_msg); return
    end
    
    shape_file_i = shape_file_i +1;
    
    for i = 1:length(shapes)
        % crop shapes
%         min_lon = min(hazard.lon);  max_lon = max(hazard.lon);
%         min_lat = min(hazard.lat);  max_lat = max(hazard.lat);        
        min_lon = min(tmp_lon);  max_lon = max(tmp_lon);
        min_lat = min(tmp_lat);  max_lat = max(tmp_lat);

        lon_ndx = shapes(i).X <= max_lon & shapes(i).X >= min_lon;
        lat_ndx = shapes(i).Y <= max_lat & shapes(i).Y >= min_lat;
        
        shapes(i).X = shapes(i).X(lon_ndx & lat_ndx);
        shapes(i).Y = shapes(i).Y(lon_ndx & lat_ndx);
        
        x_b = dx_dlon.*(shapes(i).X - lon_c);
        y_b = dy_dlat.*(shapes(i).Y - lat_c);
        if isempty(x_b) || isempty(y_b), continue; end
        try
            s(shape_file_i) = plot(x_b,y_b,'color',varargin{arg_i+2},'linewidth',varargin{arg_i+1});
        catch
            cprintf([1 0 0],err_msg); return
        end
    end
end
if ~isempty(varargin), legend(s,varargin(4:4:length(varargin)),'location','NorthWest'); legend boxoff; end


% Set axes tick labels and extent
axis([0 numel(tmp_lon) 0 numel(tmp_lat)]);

view(300,30);

set(gca,'XTick', linspace(0,numel(tmp_lon),5));
set(gca,'YTick', linspace(0,numel(tmp_lat),5));
set(gca,'XTickLabel', round(linspace(min(tmp_lon),max(tmp_lon),5).*1000)./1000);
set(gca,'YTickLabel', round(linspace(min(tmp_lat),max(tmp_lat),5).*1000)./1000);
set(gca,'ztick',[], 'zcolor', 'w');

xlabel('Longitude'); ylabel('Latitude');

set(gcf,'color', 'w');
grid off
box off

% Plot hazard
% Scale hazard x-y data
x_h = dx_dlon.*(hazard.lon - lon_c);
y_h = dy_dlat.*(hazard.lat - lat_c);

% [x_q, y_q] = meshgrid(linspace(min(x_h),max(x_h),400),linspace(min(y_h),max(y_h),400));
[x_q, y_q] = meshgrid(linspace(0,numel(tmp_lon),400),linspace(0,numel(tmp_lat),400));

hazard_Z = griddata(x_h,y_h,hazard_intensity,x_q,y_q);

hazard_Z(hazard_Z <= 0) = NaN;

% Surface for surge, contour for wind. Change to switch case statement to
% include other perils
if strcmp(hazard.peril_ID, 'TS')
    surface(x_q,y_q,hazard_Z,'FaceAlpha',0.4,'LineStyle','none','EdgeColor',[0 0 1]);
else
    freezeColors
    R = linspace(1,0.2,100);
    G = linspace(1,0.2,100);
    B = linspace(1,0.2,100);
    map = [R' G' B'];
    colormap(climada_colormap(hazard.peril_ID));
    caxis([min(hazard_intensity) max(hazard_intensity)]);
    [~,c_h] = contourf(x_q,y_q,hazard_Z);
    set(c_h,'Linestyle','none');
end
hold on

% Extras
if check_labels
    % Plot text labels for three largest losses
    [sort_vals, sort_ndx] = sort(tmp_d_a_c,'descend');
    num_lbl = 3; % specify the number of text labels
    txt_x = ones(1,num_lbl).*(0.8*numel(tmp_lon)); % the decimal factors are arbitrary label location adjustment parameters
    txt_y = ones(1,num_lbl).*(0.2*numel(tmp_lat));
    txt_v = sort_vals(1:num_lbl);
    txt_z = max(txt_v).*linspace(1.2,0.6,num_lbl);
    %Make sure factors of 10 are correct w.r.t (*)
    
    for i = 1:num_lbl
        [J, I] = find(loss_plot == sort_vals(i));
        loss_plot(sort_ndx(i)) = 0;
        lbl_y = [J(1) txt_y(i)];
        lbl_x = [I(1) txt_x(i)];
        lbl_z = [txt_v(i) txt_z(i)];
        plot3(lbl_x,lbl_y,lbl_z,'color',[0. 0.1 0.1])
    end
    text(txt_x,txt_y,txt_z,strcat(' USD ',num2str(round(txt_v./100).*100)),'fontsize',12);
end
axis tight
axis([0 numel(tmp_lon) 0 numel(tmp_lat)]);
hold off