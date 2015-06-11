function climada_distributed_measures(measure_grid, hazard)
%climada_distributed_measures

% NAME:
%   climada_distributed_measures
% PURPOSE:
%  script to account for measures that are introduced via an ascii file
%  the script reduces a selected hazard intensity by the grid value after
%  checking for NaN value and ensuring no minus values are produced. 
%
% CALLING SEQUENCE:
%   climada_distributed_measures
% EXAMPLE:
%   climada_distributed_measures
% INPUTS:
%   acsii grid
%   .mat hazard file
%   
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   new hazard file called hazard_w_measures, hazard, measure_grid (reduction from ascii
%   file)
% RESTRICTIONS:
%   none
% MODIFICATION HISTORY:
% Jacob Anz, j.anz@gmx.net, 20150609
%-mea



global climada_global
if ~climada_init_vars,return;end % init/import global variables


%open the grid measure (ascii) file
measure_grid=climada_ascii_read;


%change the ccordinates from UTM to degrees
%check if coordinates are already in UTM
if numel(num2str(abs(measure_grid.lon(1))))>8 
    
    [measure_grid.lon,measure_grid.lat]=utm2ll_shift(measure_grid.lon,measure_grid.lat);
    
  else
end


%Ask if the file is a fraction or a total value
choice = menu('Select if the ascii file contains total values (f.e. "m") or fractions ("%")','Total value','Fraction');


% prompt for hazard if not given
if ~exist('hazard','var'),hazard=[];end

if isempty(hazard) % local GUI
    hazard=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard, 'Select hazard event set (or centroids) to encode to:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard=fullfile(pathname,filename);
    end
elseif ischar(hazard)
    if strcmp(hazard,'SKIP'),return;end % special case, see climad_entity_read
end

% load the hazard, if a filename has been passed
if ~isstruct(hazard)
    hazard_file=hazard;hazard=[];
    load(hazard_file);
end


%Create the centroids
if isfield(hazard,'intensity')
    % hazard does indeed contain a hazard structure
    % hence we do not need all fields
    centroids.lon=hazard.lon;
    centroids.lat=hazard.lat;
    if isfield(hazard,'filename'),centroids.filename =hazard.filename;end
    if isfield(hazard,'comment'), centroids.comment  =hazard.comment;end
else
    % hazard does contain centroids
    centroids=hazard; clear hazard
end


% start encoding
n_points=length(centroids.lon);

t0       = clock;
msgstr   = sprintf('Encoding %i assets ... ',n_points);
mod_step = 10; % first time estimate after 10 assets, then every 100
if climada_global.waitbar
    fprintf('%s (updating waitbar with estimation of time remaining every 100th asset)\n',msgstr);
    h        = waitbar(0,msgstr);
    set(h,'Name','Encoding assets');
else
    fprintf('%s (waitbar suppressed)\n',msgstr);
    format_str='%s';
end

for asset_i=1:n_points
    if climada_global.waitbar,waitbar(asset_i/n_points,h);end
    
    dist_m=climada_geo_distance(measure_grid.lon(asset_i),measure_grid.lat(asset_i),centroids.lon,centroids.lat);
    [min_dist,min_dist_index] = min(dist_m);
   centroid_index(asset_i)=min_dist_index;
    
    %if verbose,fprintf('%f/%f --> %f/%f\n',assets.lon(asset_i),assets.lat(asset_i),centroids.lon(min_dist_index),centroids.lat(min_dist_index));end
    
    % the progress management
    if mod(asset_i,mod_step)==0
        mod_step          = 100;
        t_elapsed_event   = etime(clock,t0)/asset_i;
        events_remaining  = n_points-asset_i;
        t_projected_sec   = t_elapsed_event*events_remaining;
        if t_projected_sec<60
            msgstr = sprintf('est. %3.0f sec left (%i/%i assets)',t_projected_sec,asset_i,n_points);
        else
            msgstr = sprintf('est. %3.1f min left (%i/%i assets)',t_projected_sec/60,asset_i,n_points);
        end
        if climada_global.waitbar
            waitbar(asset_i/n_points,h,msgstr); % update waitbar
        else
            fprintf(format_str,msgstr); % write progress to stdout
            format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
        end
    end
    
end % asset_i
if climada_global.waitbar
    close(h) % dispose waitbar
else
    fprintf(format_str,''); % move carriage to begin of line
end


%Check if the the reduction contains 'ans', if it does, turn it to 0 or in
%case of a fraction turn it to 1

measure_grid.value_2=measure_grid.value;

%check if the number and coordinates of centroids match
hazard_w_measures=hazard;

%calculte the new intensity, assuming the measure redution values are
%negative
if choice ==1; 
    measure_grid.value_2(isnan(measure_grid.value_2))=0;
for i=1:length(hazard.frequency)
    
hazard_w_measures.intensity(i,:)=full(hazard.intensity(i,:))+ measure_grid.value_2;
end

else 
    measure_grid.value_2(isnan(measure_grid.value_2))=1;
for i=1:length(hazard.frequency)
    for n=1:length(hazard.intensity(1,:))
hazard_w_measures.intensity(i,n)=full(hazard.intensity(i,n))* measure_grid.value_2(n);    
    end
end
end

hazard_w_measures.measurestatement= measure_grid.comment;  

%check for negative intensities and set them to 0
     
       [find1,find2]=find(full(hazard_w_measures.intensity<0));
        hazard_w_measures.intensity(find1,find2)=0;

%plot the difference of the hazard sets without and with measures

hazard_diff=hazard;
hazard_diff.intensity=hazard.intensity-hazard_w_measures.intensity;
figure(1)
climada_hazard_plot(hazard_diff);
title('Area affected by the measures')
figure(2)
scatter3(measure_grid.lon,measure_grid.lat,measure_grid.value)
disp('done')

%saving to workspace
assignin('base','hazard_w_measures',hazard_w_measures)
assignin('base','measure_grid',measure_grid)
assignin('base','hazard',hazard)

end
%clear the variables
%clear (ans, asset_i, centroid_index, dist_m, events_remaining, filename, find1, find2,h, hazard_file, i, min_dist,min_dist_index, mod_step,msgstr, n_points, pathname, t0, t_elapsed_event,t_projected_sec);