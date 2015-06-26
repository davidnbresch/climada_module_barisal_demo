function  fig = climada_waterfall_graph_2timehorizons(return_period, check_printplot, currency,varargin)
% waterfall figure, expected damage for specified return period for
% - today,
% - increase from economic growth,
% - increase from high climate change, total expected damage 2030
% for the three EDS quoted above
% NAME:
%   climada_waterfall_graph
% PURPOSE:
%   plot expected damage for specific return period
% CALLING SEQUENCE:
%   climada_waterfall_graph(EDS1, EDS2, EDS3, return_period,
%   check_printplot)
% EXAMPLE:
%   climada_waterfall_graph
% INPUTS:
%   none
% OPTIONAL INPUT PARAMETERS:
%   EDS:            three event damage sets
%                   - today
%                   - economic growth
%                   - cc combined with economic growth, future
%   return_period:  requested return period for according expected damage,or
%                   annual expted damage, prompted if not given
%   check_printplot:if set to 1, figure saved, default 0.
% OUTPUTS:
%   waterfall graph
% MODIFICATION HISTORY:
% Lea Mueller, 20110622
% Martin Heynen, 20120329
% David N. Bresch, david.bresch@gmail.com, 20130316 EDS->EDS
% David N. Bresch, david.bresch@gmail.com, 20150419 try-catch for arrow plotting
%-

global climada_global
if ~climada_init_vars, return; end

% poor man's version to check arguments
if ~exist('EDS'             ,'var'), EDS                = [];   end
if ~exist('return_period'   ,'var'), return_period      = [];   end
if ~exist('check_printplot' ,'var'), check_printplot    = 0;    end
if ~exist('currency'        ,'var'), currency           = 0;    end


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

%% parse varargin & get all expected annual damages
if ~isstruct(varargin{1})
    cprintf([1 0 0],'ERROR: first argument after ''currency'' must be an EDS struct\n')
    return
end

EDS = struct([]); scenario_names = {};damage = []; hazard_names = {}; %init
for v_i = 1:length(varargin)
    if isstruct(varargin{v_i})      
        EDS = [EDS varargin{v_i}];
        damage =  [damage, [varargin{v_i}.ED]' ];
    elseif ischar(varargin{v_i})
        scenario_names{end+1} = varargin{v_i};
    end
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

%% set default return period, 250 years
if ~exist ('return_period', 'var'), return_period = []   ; end
if isempty(return_period)
    prompt   ='Choose specific return period or annual expected damage [e.g. 1, 10, 500, AED]:';
    name     ='Return period or annual expected damage';
    defaultanswer = {'AED or 1 or 50 etc'};
    answer   = inputdlg(prompt, name, 1, defaultanswer);
    if strcmp(answer{1},'AED')
        return_period = 9999;
    else
        return_period = str2num(answer{1});
    end
elseif ischar(return_period)
    if strcmp(return_period,'AED')
        return_period = 9999;
    else
        fprintf('Please indicate return period (e.g. 1, 34, 50) or "AED"\n "%s" does not exist\n',return_period)
        return
    end
end

%init
damage = [];
value = [];
%% collect damage for all EDS for the given return period in the damage vector
for EDS_i = 1:length(EDS)
    % check if annual expected damage is requested
    if return_period == 9999
        damage(EDS_i) = EDS(EDS_i).ED;
        value (EDS_i) = EDS(EDS_i).Value;
    else
        % find index for requested return period
        r_index = EDS(EDS_i).R_fit == return_period;
        if sum(r_index)<1
            fprintf('\n--no information available for requested return period %d year -- \n',return_period)
            
            fprintf('--calculate damage for specific return period %d  --\n', return_period)
            EDS(EDS_i) = climada_EDS_stats(EDS(EDS_i), '', return_period);
            r_index    = EDS(EDS_i).R_fit == return_period;
            
        end
        damage(EDS_i) = EDS(EDS_i).damage_fit(r_index);
    end
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
damage(length(EDS)+1) = damage(length(EDS));
% +1 to add full bar for 2030
damage_count          = length(damage)+1;
% add full bar for time horizon 2030 (EDS at position 3)
damage = [damage(1:3) damage(3) damage(4:end)];
damage                = [0 damage];

% reorder values
% value        = value(index);
value(end+1) = value(end);
% add full bar for time horizon 2030 (EDS at position 3)
value        = [value(1:3) value(3) value(4:end)];
value        = [0 value];



%% figure parameters
dmg_dig = 0;
damage = damage *10^-dmg_dig;
while max(max(damage)) > 1000
    dmg_dig = dmg_dig+3;
    damage = damage/1000;
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
if isfield(EDS,'Value_total')
    TAV_nr = unique([EDS(:).Value_total]);
else
    TAV_nr = unique([EDS(:).Value]);
end

while mean(TAV_nr) > 1000
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
    case 12
        TAV_unit = 'tn';
end

% N      = -abs(floor(log10(max(TAV_nr)))-1);
%TAV_nr = round(TAV_nr*10^N)/10^N;

% fontsize_  = 8;
fontsize_  = 12;
fontsize_2 = fontsize_ - 3;
stretch    = 0.3;

%% create colormap
% yellow - red color scheme
cmap = climada_colormap('waterfall', length(EDS));


%% create figure
% fig        = climada_figuresize(0.57,0.7);
fig        = climada_figuresize(0.57,0.9);
hold on
area([damage_count-stretch damage_count+stretch], damage(end-1)*ones(1,2),'facecolor',cmap(end-1,:),'edgecolor','none')
for i = 1:length(damage)-2
    indx = i;
    h(i) = patch( [i-stretch i+stretch i+stretch i-stretch],...
          [damage(indx) damage(indx) damage(i+1) damage(i+1)],...
          cmap(i,:),'edgecolor','none');
end
% full bar for 2030 time horizon
i = 4;
h(i) = patch( [i-stretch i+stretch i+stretch i-stretch],...
          [0 0 damage(i+1) damage(i+1)],...
          cmap(i,:),'edgecolor','none');
      
for i = 1:length(damage)-2
    if i <=3 %first time horizon
        plot([i+stretch 4-stretch],[damage(i+1) damage(i+1)],':','color',cmap(end,:))
    else % second time horizon
        plot([i+stretch damage_count-stretch],[damage(i+1) damage(i+1)],':','color',cmap(end,:))
    end
end

%% display damage string above bars
damage_disp      = diff(damage);
damage_disp(end) = damage(end);

%number of digits before the comma (>10) or behind the comma (<10)
N = 2;

%display damage string above bars
strfmt = ['%2.' int2str(N) 'f'];
dED    = 0.0;
for d_i = 2:damage_count-1
    if d_i>2 && value(d_i)<value(d_i+1)
       indx = find(value(d_i)>value);
       indx = indx(end)+1;
    else
        indx = d_i;
    end
    if d_i ~= 4
        textstr = num2str(damage(d_i+1)-damage(indx),strfmt);
        text(d_i-dED, damage(indx)+ (damage(d_i+1)-damage(indx))/2, ...
            textstr, 'color','w', 'HorizontalAlignment','center',...
            'VerticalAlignment','middle','FontWeight','bold','fontsize',fontsize_);
    end
end
text(1, damage(2)             , num2str(damage_disp(1)  ,strfmt), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);
text(damage_count, damage(end), num2str(damage_disp(end),strfmt), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);
text(4, damage(4), num2str(damage(4),strfmt), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);


%remove xlabels and ticks
set(gca,'xticklabel',[],'FontSize',10,'XTick',zeros(1,0),'layer','top');

%axis range and ylabel
xlim([0.5 damage_count+1-0.5])
ylim([0   max(damage)*1.25])
ylabel(['Damage amount (' currency ' ' dmg_unit ')'],'fontsize',fontsize_+2)



%% display arrows
% % dED2 = 0.05;% dED3 = 0.10;
% dED2 = stretch+0.05;
% dED3 = stretch+0.07;
% for d_i=2:damage_count-1
%     try
%         if d_i>2 & value(d_i)<value(d_i+1)
%            indx = find(value(d_i)>value);
%            indx = indx(end)+1;
%         else
%             indx = d_i;
%         end
%         climada_arrow ([d_i+dED2 damage(indx)], [d_i+dED2 damage(d_i+1)],...
%                        40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[0.5 0.5 0.5]);
%         %climada_arrow ([d_i+dED2 damage(d_i)], [d_i+dED2 damage(d_i+1)],...
%         %               40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[0.5 0.5 0.5]);
%         textstr = ['+' int2str((damage(d_i+1)-damage(indx))/damage(indx)*100) '%'];
%     catch
%         fprintf('Warning: arrow printing failed in %s (1)\n',mfilename);
%     end
%     text(d_i+dED3, damage(indx)+(damage(d_i+1)-damage(indx))*0.5,textstr, ...
%           'color',[0. 0. 0.],'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize_-1);
% end
% %arrow for total damage
% try
%     climada_arrow ([damage_count damage(2)], [damage_count damage(end)],...
%         40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[256 256 256]/256);
% catch
%     fprintf('Warning: arrow printing failed in %s (3)\n',mfilename);
% end
% text(damage_count, damage(2)-max(damage)*0.02, ['+' int2str((damage(end)-damage(2))/damage(2)*100) '%'],...
%     'color','w','HorizontalAlignment','center','VerticalAlignment','top','fontsize',fontsize_);


% %arrow cc
% try
%     climada_arrow ([3+dED2 damage(3)], [3+dED2 damage(4)], 40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[0.5 0.5 0.5]);
% catch
%     fprintf('Warning: arrow printing failed in %s (2)\n',mfilename);
% end
% text (3+dED3, damage(3)+diff(damage(3:4))*0.5, ['+' int2str((damage(4)-damage(3))/damage(2)*100) '%'], ...
%     'color',[0. 0. 0.],'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize_-1);
% %arrow total
% try
%     climada_arrow ([4 damage(2)], [4 damage(4)], 40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[256 256 256]/256);
% catch
%     fprintf('Warning: arrow printing failed in %s (3)\n',mfilename);
% end
% text (4, damage(2)-max(damage)*0.02, ['+' int2str((damage(4)-damage(2))/damage(2)*100) '%'], 'color','w','HorizontalAlignment','center','VerticalAlignment','top','fontsize',fontsize_);




%% set title
if return_period == 9999
    textstr = 'Annual Expected Damage (AED)';
else
    textstr = ['Expected damage with a return period of ' int2str(return_period) ' years'];
end
textstr_TAV_present = sprintf('Total asset value (%s): %s %2.0f %s',...
    num2str(min([EDS.reference_year])),currency,min(TAV_nr),TAV_unit);
textstr_TAV_future = sprintf('Total asset value (%s): %s %2.0f %s',...
    num2str(median(unique([EDS.reference_year]))),currency,median(TAV_nr),TAV_unit);
textstr_TAV_future2 = sprintf('Total asset value (%s): %s %2.0f %s',...
    num2str(max([EDS.reference_year])),currency,max(TAV_nr),TAV_unit);
if strcmpi(currency,'PEOPLE')
    textstr = 'Annual Expected no. of Casualties';
    textstr_TAV_present = sprintf('Total population (%s): %2.0f %s %s',...
        num2str(min([EDS.reference_year])),min(TAV_nr),TAV_unit,currency);
    textstr_TAV_future = sprintf('Total population (%s): %2.0f %s %s',...
        num2str(median(unique([EDS.reference_year]))),median(TAV_nr),TAV_unit,currency);
    textstr_TAV_future2 = sprintf('Total population (%s): %2.0f %s %s',...
        num2str(max([EDS.reference_year])),max(TAV_nr),TAV_unit,currency); 
end

text(1-stretch, max(max(damage))*1.21,textstr, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold','fontsize',fontsize_);
text(1-stretch, max(max(damage))*1.15,textstr_TAV_present, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','normal','fontsize',fontsize_2);
text(1-stretch, max(max(damage))*1.11,textstr_TAV_future, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','normal','fontsize',fontsize_2);
text(1-stretch, max(max(damage))*1.07,textstr_TAV_future2, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','normal','fontsize',fontsize_2);

clear textstr

%% set xlabel 
for d_i = 1:damage_count
    switch d_i
        case 1
            textstr = {[num2str(climada_global.present_reference_year) ' today''s'];'expected damage'};
        case 2
            textstr = {'Increase'; 'from econ.'; 'growth'};
        case 3
            textstr = {'Increase'; 'from'; 'climate'; 'change'};
        case 4
            textstr = {'Total';'expected';sprintf('damage %d',EDS(d_i-1).reference_year)};
        case 5
            textstr = {'Increase'; 'from econ.'; 'growth'; sprintf('%d',EDS(d_i).reference_year)};
        case 6
            textstr = {'Increase'; 'from'; 'climate'; 'change'};
        case 7
            textstr = {'Total';'expected';sprintf('damage %d',EDS(d_i-2).reference_year)};
    end
    text(d_i-stretch, damage(1)-max(damage)*0.02, textstr,...
         'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
end

%% set xlabel old
% for d_i = 2:damage_count-1
%     % economic growth (same hazard, no climate change)
%     if strcmp(EDS(d_i).hazard.filename,EDS(d_i-1).hazard.filename) & EDS(d_i).Value ~= EDS(d_i-1).Value
%         textstr = {'Increase'; 'from econ.'; 'growth'; sprintf('%d',EDS(d_i).reference_year)};
%     
%     % climate change (different hazard, same asset value, no economic growth)    
%     elseif ~strcmp(EDS(d_i).hazard.filename,EDS(d_i-1).hazard.filename) & EDS(d_i).Value == EDS(d_i-1).Value
%         %textstr = {'Increase'; 'from climate'; sprintf('change; %d',EDS(d_i).reference_year)};
%         textstr = {'Increase'; 'from mod.'; 'climate change';sprintf('%d',EDS(d_i).reference_year)};
%     
%     % climate change (same hazard, same asset value, no economic growth)    
%     elseif strcmp(EDS(d_i).hazard.filename,EDS(d_i-1).hazard.filename) & EDS(d_i).Value == EDS(d_i-1).Value
%         %textstr = {'Increase'; 'from climate'; sprintf('change; %d',EDS(d_i).reference_year)};
%         textstr = {'Increase'; 'from extr.'; 'climate change';sprintf('%d',EDS(d_i).reference_year)};    
%     
%     % just any other incremental increase    
%     else
%         textstr = {'Increase'; 'until'; sprintf('%d',EDS(d_i).reference_year)};
%     end
%     text(d_i-stretch, damage(1)-max(damage)*0.02, textstr,...
%          'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% end
% % first and last xlabel
% text(1-stretch, damage(1)-max(damage)*0.02, {[num2str(climada_global.present_reference_year) ' today''s'];'expected damage'}, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% textstr = {sprintf('%d total',EDS(d_i).reference_year);'expected';'damage'};
% text(damage_count-stretch, damage(1)-max(damage)*0.02,textstr,...
%      'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);

 
% text(2-stretch, damage(1)-max(damage)*0.02, {'Incremental increase';'from economic';'gowth; 2030'},          'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% text(3-stretch, damage(1)-max(damage)*0.02, {'Incremental increase';'2050'},                                 'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% 
% % text(2-stretch, damage(1)-max(damage)*0.02, {'Incremental increase';'from economic';'gowth; no climate';'change'},          'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% % text(3-stretch, damage(1)-max(damage)*0.02, {'Incremental increase';'from climate change'},                                 'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% % text(4-stretch, damage(1)-max(damage)*0.02, {[num2str(climada_global.future_reference_year) ', total'];'expected damage'},    'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% 
% text(2-stretch, damage(1)-max(damage)*0.02, {'Incremental increase';'from economic';'gowth; 2030'},          'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% text(3-stretch, damage(1)-max(damage)*0.02, {'Incremental increase';'2050'},                                 'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
% text(4-stretch, damage(1)-max(damage)*0.02, {'2050 total';'expected damage'},    'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);




%% add legend
%L = legend(h,legend_str(index),'location','NorthOutside','fontsize',fontsize_2);
%set(L,'Box', 'off')
if damage_count<=3
    L=legend(h, legend_str(index),'Location','NorthEast');
    set(L,'Box', 'off')
    set(L,'Fontsize',fontsize_2)
end

end % climada_waterfall_graph





