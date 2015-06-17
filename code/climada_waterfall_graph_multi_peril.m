function  fig = climada_waterfall_graph_multi_peril(check_printplot,currency,varargin)
% waterfall figure, expected damage for specified return period for
% - today,
% - increase from economic growth,
% - increase from high climate change, total expected damage 2030
% for multiple hazard types (flood, cyclone, etc). Hardwired for Barisal,
% see event_selection, which lists the index of today's damage for
% different hazard types.
% NAME:
%   climada_waterfall_graph_barisal_combined
% PURPOSE:
%   plot waterfall graph based on annual expected damage for different hazard types
% CALLING SEQUENCE:
%   fig = climada_waterfall_graph_barisal_combined(EDS, check_printplot)
% EXAMPLE:
%   climada_waterfall_graph_barisal_combined
% INPUTS:
%   none
% OPTIONAL INPUT PARAMETERS:
%   EDS: multiple event damage sets that contain scenarios on today,
%   economic growth, climate change for a future time horizon, for
%   different hazard types. The vector event_selection points to today's
%   damage fordifferent hazard types.
% OUTPUTS:
%   waterfall graph
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150610, init
%-

global climada_global
if ~climada_init_vars, return; end

% poor man's version to check arguments
if ~exist('currency'       ,'var'), currency = 'USD';       end
if ~exist('check_printplot','var'), check_printplot = 0;    end
if ~exist('varargin'       ,'var'), varargin = [];          end

%% prompt for EDS if not given
if isempty(varargin) % local GUI
    EDS=[climada_global.data_dir filesep 'results' filesep '*.mat'];
    %[filename, pathname] = uigetfile(EDS, 'Select EDS:');
    [filename, pathname] = uigetfile(EDS, 'Select EDS:','MultiSelect','on');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        if iscell(filename)
            for i = 1:length(filename)
                % rename EDS to EDS1
                vars = whos('-file', fullfile(pathname,filename{i}));
                load(fullfile(pathname,filename{i}));
                %temporarily save in EDS_temp
                EDS_temp(i) = eval(vars.name);
                clear (vars.name)
            end
            EDS = EDS_temp;
        else
            EDS = fullfile(pathname,filename);
        end
    end
end

EDS = varargin{1};
for v_i = 2:length(varargin)
    EDS = [EDS varargin{v_i}];
end
%% get EDS name into legend_str
for EDS_i = 1:length(EDS)
    % identification of EDS_i
    hazard_name       = strtok(EDS(EDS_i).hazard.comment,',');
    hazard_name       = horzcat(hazard_name, ' ', int2str(EDS(EDS_i).reference_year));
    [fP, assets_name] = fileparts(EDS(EDS_i).assets.filename);
    str               = sprintf('%s | %s',assets_name, hazard_name);
    str               = strrep(str,'_',' '); % since title is LaTEX format
    str               = strrep(str,'|','\otimes'); % LaTEX format
    legend_str{EDS_i} = str;
end % EDS_i


%% reorder damage for ascending order and add last entry again
% damage = permute(damage,[1 3 2]);
% today, eco, cc, future damage
% [damage index]        = sort(damage,'ascend');


%% get all expected annual damages

damage        = [];
hazard_names  = {};
for e_i = 1:length(varargin)
    damage =  [damage, [varargin{e_i}.ED]' ];
end
for h_i = 1:length(varargin{1})
    hazard_names{end+1} = varargin{1}(h_i).hazard.comment;
end

[n_hazards,n_scenarios] = size(damage);

damage_  = damage;
damage_(:,1) = cumsum(damage(:,1));
damage_diff = cumsum(diff(damage,1,2));
damage_ = [zeros(1,n_scenarios); damage_];
for i = 2:n_scenarios
    damage_(2:end,i) = damage_diff(:,i-1)+damage_(end,i-1);
    damage_(1,i) = damage_(end,i-1);
end
damage_(2:end,end+1) = cumsum(damage(:,end));
damage_count = n_scenarios+1;

% reorder values
% value        = value(index);
% value(end+1) = value(end);
% value        = [0 value];

% climada_waterfall_graph_advanced

%% figure parameters
dmg_dig = 0;
damage_ = damage_ *10^-dmg_dig;
while max(max(damage_)) > 1000
    dmg_dig = dmg_dig+3;
    damage_ = damage_/1000;
end
switch dmg_dig
    case 3
        dmg_unit = 'thousands';
    case 6
        dmg_unit = 'millions';
    case 9
        dmg_unit = 'billions';
    case 12
        dmg_unit = 'trillions';
    otherwise
        dmg_unit = '';
end

% TAV of portfolio
TAV_dig = 0;
TAV_nr = unique([EDS(:).Value]);
while min(TAV_nr) > 1000
    TAV_dig = TAV_dig+3;
    TAV_nr = TAV_nr./1000;
end
switch TAV_dig
    case 3
        TAV_unit = 'k';
    case 6
        TAV_unit = 'm';
    case 9
        TAV_unit = 'bn';
    case 9
        TAV_unit = 'tn';
end
TAV_nr = round(unique([EDS(:).Value])*10^-TAV_dig);
N      = -abs(floor(log10(max(TAV_nr)))-1);
TAV_nr = round(TAV_nr*10^N)/10^N;

% fontsize_  = 8;
fontsize_  = 12;
fontsize_2 = fontsize_ - 3;
stretch    = 0.3;


%% create colormap
% yellow - red color scheme
% cmap = climada_colormap('waterfall', length(EDS));
cmap_tmp = climada_colormap('waterfall', n_hazards);
cmap      = [ 70 130 180; ... % steel blue
    99 184 255; ... % steel blue 1
    255 185  15; ... % darkgoldenrod1
    238 118   0; ... % darkorange
    238  64   0 ]/255; %orangered

cmap(end+1,:) = cmap_tmp(end,:); clear cmap_tmp


%% create figure
fig        = climada_figuresize(0.57,0.9);
hold on

%plot patches
h = [];
for s_i = 1:n_scenarios+1
    for h_i = 1:n_hazards
        indx = s_i;
        h(h_i) = patch( [s_i-stretch s_i+stretch s_i+stretch s_i-stretch],...
            [damage_(h_i+1,indx) damage_(h_i+1,indx) damage_(h_i,indx) damage_(h_i,indx)],...
            cmap(h_i,:),'edgecolor','none');
    end
end
% use capital letters for first letter
hazard_names_str = regexprep(hazard_names,'(\<[a-z])','${upper($1)}');
for i = 1:length(hazard_names)
    c = textscan(hazard_names_str{i},'%s');
    c = c{1};
    hazard_names_str{i} = [c{1} ' ' c{2} ' ' c{3}];
end
clear c
L = legend(h(end:-1:1),strrep(hazard_names_str(end:-1:1),'_',' '),'Position',[0.56,0.15,0.2,0.1]);
set(L,'Box', 'off')

% plot dotted lines
for s_i = 1:n_scenarios+1
    indx = s_i;
    plot([s_i+stretch damage_count-stretch],[damage_(end,s_i) damage_(end,s_i)],':','color',cmap(end,:))
end

%display damage string above bars
damage_disp      = [damage_(end,1) diff(damage_(end,:))];
damage_disp(end) = damage_(end,end);

%number of digits before the comma (>10) or behind the comma (<10)
N = 1;

%display damage string above bars
strfmt = ['%2.' int2str(N) 'f'];
dED    = 0.0;
for d_i = 2:damage_count-1
    %if d_i>2 & value(d_i)<value(d_i+1)
    %   indx = find(value(d_i)>value);
    %   indx = indx(end)+1;
    %else
    indx = d_i;
    %end
    textstr = num2str(damage_(end,d_i)-damage_(1,indx),strfmt);
    text(d_i-dED, damage_(1,d_i)+ (damage_(end,d_i)-damage_(1,indx))/2, ...
        textstr, 'color','w', 'HorizontalAlignment','center',...
        'VerticalAlignment','middle','FontWeight','bold','fontsize',fontsize_);
end
text(1,damage_(1,2), num2str(damage_disp(1)  ,strfmt), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);
text(damage_count, damage_(end,end), num2str(damage_disp(end),strfmt), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);

%remove xlabels and ticks
set(gca,'xticklabel',[],'FontSize',10,'XTick',zeros(1,0),'layer','top');

%axis range and ylabel
xlim([0.5 damage_count+1-0.5])
ylim([0   max(damage_(:))*1.25])

ylabel(['Damage ' currency ' ' dmg_unit],'fontsize',fontsize_+2)


%% display arrows
% dED2 = 0.05;% dED3 = 0.10;
dED2 = stretch+0.05;
dED3 = stretch+0.07;
for s_i = 2:n_scenarios
    try
        indx = s_i;
        climada_arrow ([s_i+dED2 damage_(end,indx-1)], [s_i+dED2 damage_(end,s_i)],...
            40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[0.5 0.5 0.5]);
        textstr = ['+' int2str((damage_(end,s_i)-damage_(end,indx-1))/damage_(end,indx-1)*100) '%'];
    catch
        cprintf([1 0.5 0],'WARNING: arrow printing failed in %s (1)\n',mfilename);
    end
    text(s_i+dED3, damage_(end,indx-1)+(damage_(end,s_i)-damage_(end,indx-1))*0.5,textstr, ...
        'color',[0. 0. 0.],'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize_-1);
end
%arrow for total damage
try
    climada_arrow ([damage_count damage_(end,1)], [damage_count damage_(end,end)],...
        40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[256 256 256]/256);
catch
    fprintf('Warning: arrow printing failed in %s (3)\n',mfilename);
end
text(damage_count, damage_(end,1)-max(damage_(:))*0.02, ...
    ['+' int2str((damage_(end,end)-damage_(end,1))/damage_(end,1)*100) '%'],...
    'color','w','HorizontalAlignment','center','VerticalAlignment','top','fontsize',fontsize_);



%% set title
textstr = 'Annual Expected Damage (AED)';
textstr_TAV_present = sprintf('Total asset value (%s): %s %d %s',...
    num2str(min([EDS.reference_year])),currency,min(TAV_nr),TAV_unit);
textstr_TAV_future = sprintf('Total asset value (%s): %s %d %s',...
    num2str(max([EDS.reference_year])),currency,max(TAV_nr),TAV_unit);
text(1-stretch, max(max(damage_))*1.21,textstr, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold','fontsize',fontsize_);
text(1-stretch, max(max(damage_))*1.15,textstr_TAV_present, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','normal','fontsize',fontsize_2);
text(1-stretch, max(max(damage_))*1.11,textstr_TAV_future, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','normal','fontsize',fontsize_2);

clear textstr
%% set xlabel
for d_i = 2:n_scenarios
    
    
    fld_i = ['fld_' num2str(d_i)];
    
    % economic growth (same hazard, no climate change)
    if d_i ==2 || ...
            strcmp(EDS(e_i).hazard.filename,EDS(e_i-1).hazard.filename) && EDS(e_i).Value ~= EDS(e_i-1).Value
        textstr.(fld_i) = {'Increase'; 'from economic'; 'growth'; sprintf('%d',climada_global.future_reference_year)};
        
        % climate change (different hazard, same asset value, no economic growth)
    elseif d_i == 3
        %~strcmp(EDS(e_i).hazard.filename,EDS(e_i-1).hazard.filename) & EDS(e_i).Value == EDS(e_i-1).Value
        %textstr = {'Increase'; 'from climate'; sprintf('change; %d',EDS(d_i).reference_year)};
        textstr.(fld_i) = {'Increase'; 'from moderate'; 'climate change';sprintf('%d',climada_global.future_reference_year)};
        
        % climate change (same hazard, same asset value, no economic growth)
    elseif strcmp(EDS(e_i).hazard.filename,EDS(e_i-1).hazard.filename) && EDS(e_i).Value == EDS(e_i-1).Value
        %textstr = {'Increase'; 'from climate'; sprintf('change; %d',EDS(d_i).reference_year)};
        textstr.(fld_i) = {'Increase'; 'from extreme'; 'climate change';sprintf('%d',EDS(e_i).reference_year)};
        
        % just any other incremental increase
    else
        textstr.(fld_i) = {'Increase'; 'until'; sprintf('%d',EDS(e_i).reference_year)};
    end
%     text(d_i-stretch, damage_(1,1)-max(damage(:))*0.02, textstr,...
%         'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
end
% first and last xlabel


textstr.fld_1 = {'Today''s';'expected damage';sprintf('%d',climada_global.present_reference_year)};

fld_total = ['fld_' num2str(n_scenarios+1)];
textstr.(fld_total) = {'Total';'expected damage';sprintf('%d',max([EDS.reference_year]))};

% text(1-stretch, damage(1,1)-max(max(damage(:,:)))*0.02, {[num2str(climada_global.present_reference_year) ' today''s'];'expected damage'}, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% text(2-stretch, damage(1,1)-max(max(damage(:,:)))*0.02, {},          'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% text(3-stretch, damage(1,1)-max(max(damage(:,:)))*0.02, {'Incremental increase';'from climate change'},                                 'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% text(damage_count-stretch, damage(1,1)-max(max(damage(1,:)))*0.02, {[num2str(climada_global.future_reference_year) ', total'];'expected damage'},    'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);

for i = 1:damage_count
    text(i-stretch+0.3, damage_(1,1)-max(damage_(:))*0.02,textstr.(['fld_' num2str(i)]),...
        'color','k','HorizontalAlignment','center','VerticalAlignment','top','fontsize',fontsize_2);
end
% %% add legend
% L = legend(h,legend_str(index),'location','NorthOutside','fontsize',fontsize_2);
% set(L,'Box', 'off')
% if damage_count<=3
%     L=legend(h, legend_str(index),'Location','East');
%     set(L,'Box', 'off')
%     set(L,'Fontsize',fontsize_2)
% end


%% print if needed
if isempty(check_printplot)
    choice = questdlg('print?','print');
    switch choice
        case 'Yes'
            check_printplot = 1;
    end
    endw
    
    
    if check_printplot %(>=1)
        print(fig,'-dpdf',[climada_global.data_dir foldername])
        %close
        fprintf('saved 1 FIGURE in folder %s \n', foldername);
    end
    
end % climada_waterfall_graph