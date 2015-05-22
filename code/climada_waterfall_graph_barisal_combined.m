function  fig = climada_waterfall_graph_barisal_combined(EDS, check_printplot)
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
% Lea Mueller, muellele@gmail.com, 20150511, init, hard-wired for Barisal
%-

global climada_global
if ~climada_init_vars, return; end

% poor man's version to check arguments
if ~exist('EDS'            ,'var'), EDS = []; end
if ~exist('check_printplot','var'), check_printplot = 0; end


%% prompt for EDS if not given
if isempty(EDS) % local GUI
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
% load the EDS, if a filename has been passed
if ~isstruct(EDS)
    EDS_file=EDS; EDS=[];
    load(EDS_file);
    if isempty(EDS)
        if exist('EDS_all','var')
            EDS = EDS_all;
        else
            fprintf('Variable name of selected EDS does not match ''EDS''. Please check. \n')
            return
        end
    end
end


%% find different hazard types, time horizons, eco and cc scenarios
hazard_names = {'flood_depth_monsoon' 'flood_duration_monsoon' 'flood_depth_cyclone' 'flood_duration_cyclone' 'cyclone_wind'}; 

%climate change scenario
cc_scenario = {'no change' 'moderate' 'extreme'}; 

%timehorizon
% timehorizon = [2014 2030 2050];
% timehorizon = [2014 2030];
timehorizon = [2014 2050];
years       = [EDS.reference_year];
% event_selection = [1 11 6 16 21];
event_selection = [1 15 8 22 29];
    

%% get EDS name into legend_str
% for EDS_i = 1:length(EDS)
%     % identification of EDS_i
%     hazard_name       = strtok(EDS(EDS_i).hazard.comment,',');
%     hazard_name       = horzcat(hazard_name, ' ', int2str(EDS(EDS_i).reference_year));
%     [fP, assets_name] = fileparts(EDS(EDS_i).assets.filename);
%     str               = sprintf('%s | %s',assets_name, hazard_name);
%     str               = strrep(str,'_',' '); % since title is LaTEX format
%     str               = strrep(str,'|','\otimes'); % LaTEX format
%     legend_str{EDS_i} = str;
% end % EDS_i


%% reorder damage for ascending order and add last entry again
% damage = permute(damage,[1 3 2]);
% today, eco, cc, future damage
% [damage index]        = sort(damage,'ascend');


%% get all expected annual damages
ED_list       = [EDS.ED];
num_scenarios = length(cc_scenario)+length(timehorizon)-1;
damage        = zeros(length(hazard_names), num_scenarios);

for h_i = 1:length(hazard_names)
    %event_selection(h_i):event_selection(h_i)+num_scenarios-1
    damage(h_i,:) =  ED_list([event_selection(h_i) event_selection(h_i)+4 event_selection(h_i)+5 event_selection(h_i)+6]);
end

damage_prepared  = damage;
damage_prepared(:,1) = cumsum(damage(:,1));
damage_diff = cumsum(diff(damage,1,2));
damage_prepared = [zeros(1,num_scenarios); damage_prepared];
for i = 2:num_scenarios
    damage_prepared(2:end,i) = damage_diff(:,i-1)+damage_prepared(end,i-1);
    damage_prepared(1,i) = damage_prepared(end,i-1);
end
damage_prepared(2:end,end+1) = cumsum(damage(:,end));
damage_count = num_scenarios+1;

% reorder values
% value        = value(index);
% value(end+1) = value(end);
% value        = [0 value];

% climada_waterfall_graph_advanced


%% define figure parameters
%digits of damage
% digits = log10(max(damage));
% digits = floor(digits)-1;
% digits = 9;
% digits = 6;
digits = 0;
% damage = damage*10^-digits;
dig    = digits;

% TIV of portfolio
% TIV_nr = round(unique([EDS(:).Value])*10^-digits);
% N      = -abs(floor(log10(max(TIV_nr)))-1);
% TIV_nr = round(TIV_nr*10^N)/10^N;

TIV_nr = unique([EDS(event_selection).Value])*10^-digits;
% N      = -abs(floor(log10(max(TIV_nr)))-1);
% TIV_nr = round(TIV_nr*10^N)/10^N;

% fontsize_  = 8;
fontsize_  = 12;
fontsize_2 = fontsize_ - 3;
% stretch    = 0.3;
stretch    = 0.3;



%% create colormap
% yellow - red color scheme
% cmap = climada_colormap('waterfall', length(EDS));
cmap_temp = climada_colormap('waterfall', length(hazard_names));
cmap      = [ 70 130 180; ... % steel blue
              99 184 255; ... % steel blue 1
             255 185  15; ... % darkgoldenrod1
             238 118   0; ... % darkorange
             238  64   0 ]/255; %orangered
cmap(end+1,:) = cmap_temp(end,:);


%% create figure
% fig        = climada_figuresize(0.57,0.7);
fig        = climada_figuresize(0.57,0.9);
hold on

%plot patches
h = [];
for s_i = 1:num_scenarios+1
    for h_i = 2:length(hazard_names)+1
        indx = s_i;
        h(h_i) = patch( [s_i-stretch s_i+stretch s_i+stretch s_i-stretch],...
              [damage_prepared(h_i,indx) damage_prepared(h_i,indx) damage_prepared(h_i-1,indx) damage_prepared(h_i-1,indx)],...
              cmap(h_i-1,:),'edgecolor','none');
    end
end
% use capital letters for first letter
hazard_names_string = regexprep(hazard_names,'(\<[a-z])','${upper($1)}');
L = legend(h(end:-1:2),strrep(hazard_names_string(end:-1:1),'_',' '),'Position',[0.56,0.15,0.2,0.1]);
set(L,'Box', 'off')

% plot dotted lines
for s_i = 1:num_scenarios+1
    indx = s_i;
    plot([s_i+stretch damage_count-stretch],[damage_prepared(end,s_i) damage_prepared(end,s_i)],':','color',cmap(end,:))
end

%display damage string above bars
damage_disp      = [damage_prepared(end,1) diff(damage_prepared(end,:))];
damage_disp(end) = damage_prepared(end,end);

%number of digits before the comma (>10) or behind the comma (<10)
% if max(damage)>10
%     N = -abs(floor(log10(max(damage)))-1);
%     damage_disp = round(damage_disp*10^N)/10^N;
%     N = 0;
% else
    %N = round(log10(max(damage_disp)));
    N = 1;
% end


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
    textstr = num2str(damage_prepared(end,d_i)-damage_prepared(1,indx),strfmt);
    text(d_i-dED, damage_prepared(1,d_i)+ (damage_prepared(end,d_i)-damage_prepared(1,indx))/2, ...
         textstr, 'color','w', 'HorizontalAlignment','center',...
         'VerticalAlignment','middle','FontWeight','bold','fontsize',fontsize_);
end
text(1,damage_prepared(1,2), num2str(damage_disp(1)  ,strfmt), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);
text(damage_count, damage_prepared(end,end), num2str(damage_disp(end),strfmt), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);

%remove xlabels and ticks
set(gca,'xticklabel',[],'FontSize',10,'XTick',zeros(1,0),'layer','top');

%axis range and ylabel
xlim([0.5 damage_count+1-0.5])
ylim([0   max(damage_prepared(:))*1.25])
if dig == 0
    ylabel('Damage (mn BDT)','fontsize',fontsize_+2)
else
    ylabel(['Damage amount \cdot 10^{', int2str(dig) '}'],'fontsize',fontsize_+2)
end


%% display arrows
% dED2 = 0.05;% dED3 = 0.10;
dED2 = stretch+0.05;
dED3 = stretch+0.07;
for s_i = 2:num_scenarios
    try
        %if s_i>2 & value(s_i)<value(s_i+1)
        %   indx = find(value(s_i)>value);
        %   indx = indx(end)+1;
        %else
            indx = s_i;
        %end
        climada_arrow ([s_i+dED2 damage_prepared(end,indx-1)], [s_i+dED2 damage_prepared(end,s_i)],...
                       40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[0.5 0.5 0.5]);
        textstr = ['+' int2str((damage_prepared(end,s_i)-damage_prepared(end,indx-1))/damage_prepared(end,indx-1)*100) '%'];
    catch
        fprintf('Warning: arrow printing failed in %s (1)\n',mfilename);
    end
    text(s_i+dED3, damage_prepared(end,indx-1)+(damage_prepared(end,s_i)-damage_prepared(end,indx-1))*0.5,textstr, ...
          'color',[0. 0. 0.],'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize_-1);
end
%arrow for total damage
try
    climada_arrow ([damage_count damage_prepared(end,1)], [damage_count damage_prepared(end,end)],...
        40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[256 256 256]/256);
catch
    fprintf('Warning: arrow printing failed in %s (3)\n',mfilename);
end
text(damage_count, damage_prepared(end,1)-max(damage_prepared(:))*0.02, ...
    ['s+' int2str((damage_prepared(end,end)-damage_prepared(end,1))/damage_prepared(end,1)*100) '%'],...
    'color','w','HorizontalAlignment','center','VerticalAlignment','top','fontsize',fontsize_);



%% set title
textstr = 'Annual Expected Damage (AED)';
if dig == 0
    textstr_TIV_2 = sprintf('%4.0f, ', TIV_nr);
    textstr_TIV_3 = ' (mn BDT)';
else
    textstr_TIV_2 = sprintf('%d, ', TIV_nr);
    textstr_TIV_3 = sprintf('10^%d USD', digits);
end
textstr_TIV_1 = 'Total assets: ';
textstr_TIV_2(end-1:end) = [];
textstr_TIV = [textstr_TIV_1 textstr_TIV_2 textstr_TIV_3];

text(1-stretch, max(damage_prepared(:))*1.20,textstr, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold','fontsize',fontsize_);
text(1-stretch, max(damage_prepared(:))*1.15,textstr_TIV, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','normal','fontsize',fontsize_2);



%% set xlabel
e_ = [1 5 6 7];
for d_i = 2:num_scenarios
    
    %e_i = event_selection(1)+d_i-1;
    e_i = e_(d_i);
    % economic growth (same hazard, no climate change)
    if strcmp(EDS(e_i).hazard.filename,EDS(e_i-1).hazard.filename) & EDS(e_i).Value ~= EDS(e_i-1).Value
        textstr = {'Increase'; 'from economic'; 'growth'; sprintf('(%d)',EDS(e_i).reference_year)};
    
    % climate change (different hazard, same asset value, no economic growth)    
    elseif d_i == 3
        %~strcmp(EDS(e_i).hazard.filename,EDS(e_i-1).hazard.filename) & EDS(e_i).Value == EDS(e_i-1).Value
        %textstr = {'Increase'; 'from climate'; sprintf('change; %d',EDS(d_i).reference_year)};
        textstr = {'Increase'; 'from moderate'; 'climate change';sprintf('(%d)',EDS(e_i).reference_year)};
    
    % climate change (same hazard, same asset value, no economic growth)    
    elseif strcmp(EDS(e_i).hazard.filename,EDS(e_i-1).hazard.filename) & EDS(e_i).Value == EDS(e_i-1).Value
        %textstr = {'Increase'; 'from climate'; sprintf('change; %d',EDS(d_i).reference_year)};
        textstr = {'Increase'; 'from extreme'; 'climate change';sprintf('(%d)',EDS(e_i).reference_year)};    
    
    % just any other incremental increase    
    else
        textstr = {'Increase'; 'until'; sprintf('%d',EDS(e_i).reference_year)};
    end
    text(d_i-stretch, damage_prepared(1,1)-max(damage(:))*0.02, textstr,...
         'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
end
% first and last xlabel
textstr = {'Today''s';'expected damage';sprintf('(%d)',climada_global.present_reference_year)};
text(1-stretch, damage_prepared(1,1)-max(damage_prepared(:))*0.02, textstr, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
textstr = {'Total';'expected damage';sprintf('(%d)',EDS(e_(end)).reference_year)};
text(damage_count-stretch, damage_prepared(1,1)-max(damage_prepared(:))*0.02,textstr,...
     'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);


%% add legend
%L = legend(h,legend_str(index),'location','NorthOutside','fontsize',fontsize_2);
%set(L,'Box', 'off')
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
end


if check_printplot %(>=1)
    print(fig,'-dpdf',[climada_global.data_dir foldername])
    %close
    fprintf('saved 1 FIGURE in folder %s \n', foldername);
end

end % climada_waterfall_graph





