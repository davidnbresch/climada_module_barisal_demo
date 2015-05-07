function  fig = climada_waterfall_graph_barisal(EDS, return_period, check_printplot)
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
% Lea MuEDler, 20110622
% Martin Heynen, 20120329
% David N. Bresch, david.bresch@gmail.com, 20130316 EDS->EDS
% David N. Bresch, david.bresch@gmail.com, 20150419 try-catch for arrow plotting
%-

global climada_global
if ~climada_init_vars, return; end

% poor man's version to check arguments
if ~exist('EDS'           ,'var'), EDS = []; end
% if ~exist('EDS2'           ,'var'), EDS2 = []; end
% if ~exist('EDS3'           ,'var'), EDS3 = []; end
if ~exist('EDS_comparison','var'),EDS_comparison='';end
if ~exist('return_period'  ,'var'), return_period   = []; end
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
end



%% set default return period, 250 years
if ~exist ('return_period', 'var'), return_period = []   ; end
%if isempty(return_period)         , return_period = 9999 ; end
%if isempty(return_period)         , return_period = 10 ; end
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
            
            %R_show = sprintf('%d, ', EDS(EDS_i).R_fit');
            %R_show(end-1:end) = [];
            %fprintf('--please sEDect one of the following return periods:  --\n %s\n', R_show)
            %damage = [];
            %return
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
damage_count          = length(damage);
damage                = [0 damage];

% reorder values
% value        = value(index);
value(end+1) = value(end);
value        = [0 value];



%% define figure parameters
%digits of damage
% digits = log10(max(damage));
% digits = floor(digits)-1;
% digits = 9;
% digits = 6;
digits = 0;
damage = damage*10^-digits;
dig    = digits;

% TIV of portfolio
% TIV_nr = round(unique([EDS(:).Value])*10^-digits);
% N      = -abs(floor(log10(max(TIV_nr)))-1);
% TIV_nr = round(TIV_nr*10^N)/10^N;

TIV_nr = unique([EDS(:).Value])*10^-digits;
% N      = -abs(floor(log10(max(TIV_nr)))-1);
% TIV_nr = round(TIV_nr*10^N)/10^N;

% fontsize_  = 8;
fontsize_  = 12;
fontsize_2 = fontsize_ - 3;
% stretch    = 0.3;
stretch    = 0.3;



%% create colormap
% yellow - red color scheme
cmap = climada_colormap('waterfall', length(EDS));
% if not only today, eco, climate change, but a multiple of EDS, create 
% colormap accoringly and add last color for grey dotted line in between
% if length(EDS)>3 
%     color_line = cmap(end,:);
%     cmap = jet(length(EDS));
%     cmap = [cmap; color_line];
% end

% % green color scheme
% color_     = [227 236 208;...   %today
%               194 214 154;...   %eco
%               181 205  133;...  %clim
%               197 190 151;...   %total risk
%               120 120 120]/256; %dotted line]/255;
% color_(1:4,:) = brighten(color_(1:4,:),-0.5);



%% create figure
% fig        = climada_figuresize(0.57,0.7);
fig        = climada_figuresize(0.57,0.9);
hold on
area([damage_count-stretch damage_count+stretch], damage(end-1)*ones(1,2),'facecolor',cmap(end-1,:),'edgecolor','none')
for i = 1:length(damage)-2
    %if i>2 & value(i)<value(i+1)
    %    indx = find(value(i)>value);
    %    indx = indx(end)+1;
    %else
        indx = i;
    %end
    h(i) = patch( [i-stretch i+stretch i+stretch i-stretch],...
          [damage(indx) damage(indx) damage(i+1) damage(i+1)],...
          cmap(i,:),'edgecolor','none');
end
for i = 1:length(damage)-2
    if i==1
        plot([i+stretch damage_count+stretch],[damage(i+1) damage(i+1)],':','color',cmap(end,:))
    else
%         if i>2 & value(i)<value(i+1)
%             
%         else
%             lkjd
%         end
        plot([i+stretch damage_count-stretch],[damage(i+1) damage(i+1)],':','color',cmap(end,:))
    end
end

%display damage string above bars
damage_disp      = diff(damage);
damage_disp(end) = damage(end);

%number of digits before the comma (>10) or behind the comma (<10)
% if max(damage)>10
%     N = -abs(floor(log10(max(damage)))-1);
%     damage_disp = round(damage_disp*10^N)/10^N;
%     N = 0;
% else
    %N = round(log10(max(damage_disp)));
    N = 2;
% end


%display damage string above bars
strfmt = ['%2.' int2str(N) 'f'];
dED    = 0.0;
for d_i = 2:damage_count-1
    text(d_i-dED, damage(d_i)+ (damage(d_i+1)-damage(d_i))/2, ...
         num2str(damage_disp(d_i),strfmt), 'color','w', 'HorizontalAlignment','center',...
         'VerticalAlignment','middle','FontWeight','bold','fontsize',fontsize_);
end
text(1, damage(2)             , num2str(damage_disp(1)  ,strfmt), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);
text(damage_count, damage(end), num2str(damage_disp(end),strfmt), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);

% %damages above barsn -- int2str
% dED = 0.0;
% text(1, damage(2)                     , int2str(damage_disp(1)), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);
% text(2-dED, damage(2)+ (damage(3)-damage(2))/2, int2str(damage_disp(2)), 'color','w', 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontWeight','bold','fontsize',fontsize_);
% text(3-dED, damage(3)+ (damage(4)-damage(3))/2, int2str(damage_disp(3)), 'color','w', 'HorizontalAlignment','center', 'VerticalAlignment','middle','FontWeight','bold','fontsize',fontsize_);
% text(4, damage(4)                     , int2str(damage_disp(4)), 'color','k', 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontWeight','bold','fontsize',fontsize_);

%remove xlabels and ticks
set(gca,'xticklabel',[],'FontSize',10,'XTick',zeros(1,0),'layer','top');

%axis range and ylabel
xlim([0.5 damage_count+1-0.5])
ylim([0   max(damage)*1.25])
if dig == 0
    ylabel('Damage (1000 BDT)','fontsize',fontsize_+2)
else
    ylabel(['Damage amount \cdot 10^{', int2str(dig) '}'],'fontsize',fontsize_+2)
end


%% display arrows
% dED2 = 0.05;% dED3 = 0.10;
dED2 = stretch+0.05;
dED3 = stretch+0.07;
for d_i=2:damage_count-1
    try
        climada_arrow ([d_i+dED2 damage(d_i)], [d_i+dED2 damage(d_i+1)],...
                       40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[0.5 0.5 0.5]);
    catch
        fprintf('Warning: arrow printing failed in %s (1)\n',mfilename);
    end
    text (d_i+dED3, damage(d_i)+diff(damage(d_i:d_i+1))*0.5, ['+' int2str((damage(d_i+1)-damage(d_i))/damage(d_i)*100) '%'], ...
          'color',[0. 0. 0.],'HorizontalAlignment','left','VerticalAlignment','middle','fontsize',fontsize_-1);
end
%arrow for total damage
try
    climada_arrow ([damage_count damage(2)], [damage_count damage(end)],...
        40, 10, 30,'width',1.5,'Length',10, 'BaseAngle',90, 'EdgeColor','none', 'FaceColor',[256 256 256]/256);
catch
    fprintf('Warning: arrow printing failed in %s (3)\n',mfilename);
end
text(damage_count, damage(2)-max(damage)*0.02, ['+' int2str((damage(end)-damage(2))/damage(2)*100) '%'],...
    'color','w','HorizontalAlignment','center','VerticalAlignment','top','fontsize',fontsize_);


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
if dig == 0
    textstr_TIV_2 = sprintf('%4.0f, ', TIV_nr);
    textstr_TIV_3 = ' (1000 BDT)';
else
    textstr_TIV_2 = sprintf('%d, ', TIV_nr);
    textstr_TIV_3 = sprintf('10^%d USD', digits);
end
textstr_TIV_1 = 'Total assets: ';
textstr_TIV_2(end-1:end) = [];
textstr_TIV = [textstr_TIV_1 textstr_TIV_2 textstr_TIV_3];

text(1-stretch, max(damage)*1.20,textstr, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','bold','fontsize',fontsize_);
text(1-stretch, max(damage)*1.15,textstr_TIV, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','FontWeight','normal','fontsize',fontsize_2);



%% set xlabel
for d_i = 2:damage_count-1
    % economic growth (same hazard, no climate change)
    if strcmp(EDS(d_i).hazard.filename,EDS(d_i-1).hazard.filename) & EDS(d_i).Value ~= EDS(d_i-1).Value
        textstr = {'Increase'; 'from econ.'; 'growth'; sprintf('%d',EDS(d_i).reference_year)};
    
    % climate change (different hazard, same asset value, no economic growth)    
    elseif ~strcmp(EDS(d_i).hazard.filename,EDS(d_i-1).hazard.filename) & EDS(d_i).Value == EDS(d_i-1).Value
        %textstr = {'Increase'; 'from climate'; sprintf('change; %d',EDS(d_i).reference_year)};
        textstr = {'Increase'; 'from mod.'; 'climate change';sprintf('%d',EDS(d_i).reference_year)};
    
    % climate change (same hazard, same asset value, no economic growth)    
    elseif strcmp(EDS(d_i).hazard.filename,EDS(d_i-1).hazard.filename) & EDS(d_i).Value == EDS(d_i-1).Value
        %textstr = {'Increase'; 'from climate'; sprintf('change; %d',EDS(d_i).reference_year)};
        textstr = {'Increase'; 'from extr.'; 'climate change';sprintf('%d',EDS(d_i).reference_year)};    
    
    % just any other incremental increase    
    else
        textstr = {'Increase'; 'until'; sprintf('%d',EDS(d_i).reference_year)};
    end
    text(d_i-stretch, damage(1)-max(damage)*0.02, textstr,...
         'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
end
% first and last xlabel
text(1-stretch, damage(1)-max(damage)*0.02, {[num2str(climada_global.present_reference_year) ' today''s'];'expected damage'}, 'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);
textstr = {sprintf('%d total',EDS(d_i).reference_year);'expected';'damage'};
text(damage_count-stretch, damage(1)-max(damage)*0.02,textstr,...
     'color','k','HorizontalAlignment','left','VerticalAlignment','top','fontsize',fontsize_2);

 
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





