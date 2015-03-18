function hazard = climada_ma_hazard_set(years, centroids, hazard_set_file, check_plots)
% climada MA hazard event set
% NAME:
%   climada_ma_hazard_set
% PURPOSE:
%   Construct hazard set structure from historical daily precipitation data
%   for Asian Monsoon - specifically from the APHRODITE data, see read_APHRO_MA_V1101
% CALLING SEQUENCE:
%   hazard=climada_ma_hazard_set(years,centroids,hazard_set_file,,check_plots)
% EXAMPLE:
%   hazard=climada_ma_hazard_set
%   hazard=climada_ma_hazard_set(1983:2001,centroids,[],1)
% INPUTS:
%   hazard_set_file: the name of the newly created storm surge (TS) hazard
%       event set (if ='NO_SAVE', the hazard is just returned, not saved)
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   check_plots: whether to show plots (default = 1)
% OUTPUTS:
%   hazard: a hazard event set, see core climada doc
%       also written to a .mat file (see hazard_set_file)
%       NOTE: the monsoon asia hazard set consists only of historical events.
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150121
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('years',              'var'), years               =2007;      end
if ~exist('hazard_set_file',    'var'), hazard_set_file     ='';        end
if ~exist('centroids',          'var'), centroids           =[];        end
if ~exist('check_plots',        'var'), check_plots         =0;         end
% PARAMETERS
%
module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
if ~isdir(module_data_dir),mkdir(fileparts(module_data_dir),'data');end % create the data dir, should it not exist (no further checking)

if ~isstruct(centroids)
    centroids = climada_centroids_load;
end

% prep the region we need (rectangular region encompassing the hazard centroids)
centroids_rect=[min(centroids.lon) max(centroids.lon)...
    min(centroids.lat) max(centroids.lat)];


% prompt for hazard_set_file if not given
if isempty(hazard_set_file) % local GUI
    hazard_set_file=[module_data_dir filesep 'hazards' filesep 'test_MA_hazard.mat'];
    [fN, fP] = uiputfile(hazard_set_file, 'Save new MA hazard event set as:');
    if isequal(fN,0) || isequal(fP,0)
        return; % cancel
    else
        hazard_set_file=fullfile(fP,fN);
    end
end

precip_file = [module_data_dir filesep 'precip_data' filesep 'APHRO_MA_025deg_V1101.1951-2007.mat'];

if ~exist(precip_file,'file')
    [precip_grid, precip_array] = read_APHRO_MA_V1101(1951:2007,centroids_rect,precip_file,[],check_plots);
else
    fprintf('loading raw MA precip data from %s \n',precip_file);
    load(precip_file);
end

min_year   = min(years);
max_year   = max(years);

% select relevant time period
t_ndx = precip_array.time >= datenum(min_year,01, 01) & precip_array.time <= datenum(max_year,12, 31);

precip_array.time = precip_array.time(t_ndx);
precip_array.precip = precip_array.precip(:,t_ndx);
precip_grid.time = precip_grid.time(t_ndx);
precip_grid.precip = precip_grid.precip(:,:,t_ndx);
clear t_ndx

orig_years = max_year - min_year+1;

% fill the hazard structure
hazard.reference_year   = climada_global.present_reference_year; % default for present hazard is normally 2010
hazard.lon              = centroids.lon;
hazard.lat              = centroids.lat;
hazard.centroid_ID      = centroids.centroid_ID;
if isfield(centroids,'elevation_m'),hazard.elevation_m=centroids.elevation_m;end
hazard.orig_years       = orig_years;
hazard.orig_event_count = length(precip_grid.time);
hazard.event_count      = length(precip_grid.time);
hazard.event_ID         = 1:hazard.event_count;
hazard.orig_event_flag  = ones(1,hazard.event_count);
hazard.yyyy             = str2num(datestr(precip_grid.time,'yyyy'));
hazard.mm               = str2num(datestr(precip_grid.time,'mm'));
hazard.dd               = str2num(datestr(precip_grid.time,'dd'));
hazard.datenum          = precip_grid.time;

event_frequency         = 1/(orig_years);
hazard.frequency        = ones(1,hazard.event_count)*event_frequency;

hazard.peril_ID='MA';
hazard.comment=sprintf('MA hazard event set, generated %s',datestr(now));

try
    hazard.intensity        = spalloc(hazard.event_count,length(hazard.lon),ceil(hazard.event_count*length(hazard.lon)*0.3));
catch
    hazard.intensity        = zeros(hazard.event_count,numel(centroids.centroid_ID));
end
fprintf('processing MA precipitation at centroids for %i events...\n',hazard.event_count)
mod_step = 10; format_str = '%s'; t0 = clock;
% [LON, LAT] = meshgrid(double(precip_grid.lon),double(precip_grid.lat));
LON = double(precip_grid.lon);
LAT = double(precip_grid.lat);
for event_i = 1:hazard.event_count
    if any(any(precip_grid.precip(:,:,event_i)))
        % hazard.intensity(event_i,arr_i)=qinterp2(LON, LAT,squeeze(double(precip_grid.precip(:,:,event_i))'),hazard.lon,hazard.lat);
        % hazard.intensity(event_i,arr_i)=interp2(LON, LAT, squeeze(double(precip_grid.precip(:,:,event_i))'),hazard.lon(arr_i),hazard.lat(arr_i));
        hazard.intensity(event_i,:)=interp2(LON, LAT, squeeze(double(precip_grid.precip(:,:,event_i))'),hazard.lon(:),hazard.lat(:));
    end
    % the progress management
    if mod(event_i,mod_step)==0
        mod_step          = 100;
        t_elapsed_event   = etime(clock,t0)/event_i;
        events_remaining  = hazard.event_count-event_i;
        t_projected_sec   = t_elapsed_event*events_remaining;
        if t_projected_sec<60
            msgstr = sprintf('est. %3.0f sec left (%i/%i events)',t_projected_sec,event_i,hazard.event_count);
        else
            msgstr = sprintf('est. %3.1f min left (%i/%i events)',t_projected_sec/60,event_i,hazard.event_count);
        end
        fprintf(format_str,msgstr); % write progress to stdout
        format_str=[repmat('\b',1,length(msgstr)) '%s']; % back to begin of line
    end
end
try 
    hazard.intensity    = sparse(hazard.intensity);
end
fprintf(format_str,sprintf('processing rainfall at %i centroids for %i events took %3.1f seconds \n', ...
    numel(hazard.centroid_ID),hazard.event_count,etime(clock,t0)));

if isfield(hazard,'filename'),hazard.filename_source=hazard.filename;end
hazard.filename=hazard_set_file;
hazard.date=datestr(now);
hazard.matrix_density=nnz(hazard.intensity)/numel(hazard.intensity);
hazard.units='mm'; % store the SI unit of the hazard intensity
if ~isfield(hazard,'orig_event_count') % fix a minor issue with some hazard sets
    if isfield(hazard,'orig_event_flag')
        fprintf('field hazard.orig_event_count inferred from hazard.orig_event_flag\n')
        hazard.orig_event_count=sum(hazard.orig_event_flag);
    else
        fprintf('WARNING: no field hazard.orig_event_flag\n')
    end
end

if ~strcmp(hazard_set_file,'NO_SAVE');
    fprintf('saving MA monsoon Asia hazard set as %s\n',hazard_set_file);
    try
        save(hazard_set_file,'hazard','-v7.3');
    catch
        cprintf([1 0 0], 'ERROR: can not write to file, try saving MA hazard manually \n')
    end
end

if check_plots,figure('color','w'); climada_hazard_plot_hr(hazard,0);end % show max rainfall over ALL events

return
