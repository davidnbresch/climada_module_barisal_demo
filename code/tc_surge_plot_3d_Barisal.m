function tc_surge_plot_3d_Barisal(hazard,event_i,seashore,EDS)
% climada
% NAME:
%   tc_surge_plot_3d_Barisal
% PURPOSE:
%   plot surge hazard as 3D surface, 'overshadowing' terrain
% CALLING SEQUENCE:
%   tc_surge_plot3d(hazard,event_i);
% EXAMPLE:
%   climada_hazard_plot(climada_hazard_load);
% INPUTS:
%   hazard: hazard structure
%       consider to use climada_hazard_load, as in the example
%   event_i: the i-th event in the hazard event set to be displayed
%       if event_i=0, the maximum intensity at each centroid is shown
%       if event_i=-i, the i-th 'largest' event (sum of intensities) is shown
%           e.g. for event_i=-2, the second largest event is shown
%       default=-1 (just to get something on the screen ;-)
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   figure
% MODIFICATION HISTORY:
% David N. Bresch, david.bresch@gmail.com, 20140503
% Gilles Stassen, gillesstassen@hotmail.com 20141126 adapted
% tc_surge_plot_3d with some Barisal specific modifications.
%-

%global climada_global
%if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('hazard','var'),return;end
if ~exist('event_i','var'),event_i=-1;end
if ~exist('seashore','var'),seashore = [];end
% PARAMETERS
%
% whether we show the seashore (=1) or not (=0)
show_seashore=1; % default=1
%
% whether we show the mean sea level (MSL) as surface (=1) or not (=0)
show_MSL=0; % rather not, =0
%
% mask deep(er) sea and high(er) elevations (not of interest here) (=1)
mask_elev=1; % default=1
%
% plot centroids (to check for 'mock resolution).
plot_centroids=0; %=1 (yes) or =0 (no)


if event_i<0
    % search for i-thlargest event
    event_sum=sum(hazard.intensity,2);
    [~,sorted_i]=sort(event_sum);
    event_ii=sorted_i(length(sorted_i)+event_i+1);
    surge_height=full(hazard.intensity(event_ii,:)); % extract one event
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
    surge_height=full(max(hazard.intensity)); % max intensity at each point
    title_str=sprintf('%s max intensity at each centroid',hazard.peril_ID);
else
    surge_height=full(hazard.intensity(event_i,:)); % extract one event
    title_str=sprintf('%s event %i',hazard.peril_ID,event_i);
    % plot some further info to sdout:
    if (isfield(hazard,'name') && isfield(hazard,'yyyy')) && (isfield(hazard,'mm') && isfield(hazard,'dd'))
        fprintf('%s, %4.4i%2.2i%2.2i, event %i\n',hazard.name{event_i},hazard.yyyy(event_i),hazard.mm(event_i),hazard.dd(event_i),event_i);
    end
end
if isfield(hazard,'units'),title_str=[title_str ' [' hazard.units ']'];end % add units

if isfield(hazard,'elev')
    surge_height=surge_height+hazard.elev; % add elevation, since surge height is relative to terrain
else
    fprintf('WARNING: no elevation at centroids found, surge surface might not be correct\n');
end

% define the area of interest
bathy_coords=[min(hazard.lon) max(hazard.lon) min(hazard.lat) max(hazard.lat)];

% NOTE: here, one might set a smaller/larger area, it usually works
%bathy_coords=[90 91 22.5 23.5]

% get the bathymetry
BATI=etopo_get(bathy_coords);
X=BATI.x; % just copied for convenience matters (yes, uses memory)
Y=BATI.y;
Z=BATI.h;

% mask deep(er) sea and high(er) elevations (not of interest here)
if mask_elev
    pos_excess_height = Z >25;
    Z(pos_excess_height)=NaN;
end % mask_elev

figure('Name',title_str,'Color',[1 1 1]);
fprintf('%s\n',title_str);

% plot bathymetry
% ---------------
map = [
    40  54 154
    0 201  50
    30 211 104
    94 224 116
    162 235 130
    223 248 146
    246 229 149
    200 178 118
    162 126  94
    143  97  84
    162 125 116
    178 150 139
    199 176 170
    219 205 202
    236 228 226
    255 255 255
    ]./255;
colormap(map);

surface(X,Y,Z,'FaceAlpha',1,'LineStyle','none');
hold on;
colorbar;
caxis([-1 15]); % a reasonable colorbar for elevation

if show_MSL
    surf=surface(X,Y,Z*0,'FaceAlpha',0.1,'LineStyle','none','EdgeColor',[0 0 1]);
    % change color
    CData=get(surf,'CData');
    CData=CData*0-5; % set to a reasonable color value
    set(surf,'CData',CData);
end % show_MSL

if show_seashore
    if ~isempty(seashore)
        plot(seashore(:,1),seashore(:,2),'color',[81 81 81]/256,'linewidth',2)
    else
        climada_plot_world_borders(1,[],[],1)
    end
    axis(bathy_coords);
end % show_seashore

% overlay storm surge field
% -------------------------

% grid back to display resolution
surge_Z = griddata(hazard.lon,hazard.lat,surge_height,X,Y);
if mask_elev,surge_Z(pos_excess_height) = NaN; end

surf=surface(X,Y,surge_Z,'FaceAlpha',0.4,'LineStyle','none','EdgeColor',[0 0 1]);
% change color
CData=get(surf,'CData');
CData=CData*0-20; % set to a reasonable color value
set(surf,'CData',CData);
axis([bathy_coords(1) bathy_coords(2) bathy_coords(3) bathy_coords(4)]);

if plot_centroids
    plot3(hazard.lon,hazard.lat,hazard.lat+5,'.r','MarkerSize',1);
    fprintf('very tiny red dots: hazard centroids\n');
end

% Attempt at plotting houses....
if exist('EDS','var')
    hx = EDS.assets.lon(EDS.assets.Value>0);
    hy = EDS.assets.lat(EDS.assets.Value>0);
    
    asset_elev = griddata(hazard.lon,hazard.lat,hazard.elev,hx,hy);
    hz = asset_elev;
   
    scatter3(hx,hy,hz,'s','filled','k')
    % create 'houses'
%     vertices = [  -0.5 -0.5 -0.5; -0.5 0.5 -0.5; -0.5 0.5 0.5; -0.5 -0.5 0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5;0.5 0.5 0.5; 0.5 -0.5 0.5];
%     faces = [ 1 2 3 4;5 6 7 8; 2 6 7 3; 1 5 8 4; 1 2 6 5; 4 3 7 8];
%     
%     vertices = 0.1*vertices;
%     for i = 1 : numel(asset_elev)
%         vertices = vertices + repmat([hx(i) hy(i) hz(i)+5],size(vertices,1),1);
%         h = patch('faces',faces,'vertices',vertices);
%         hold on
%         set(h,'Facecolor',[1 0 0])
%     end
     
end
set(gca,'ztick',[]);
view(300,30);

dx=max(BATI.x(2,2)-BATI.x(2,1),BATI.x(2,2)-BATI.x(1,2)); % in degree
dx=40000/360*cos(bathy_coords(1)/90*pi/2)*dx; % in km, hard-wired for ETOPO1
fprintf('voxel extent 1'', ~%2.1f km (area ~%2.1f km^2)',dx,dx^2);
fprintf('\n');

return
