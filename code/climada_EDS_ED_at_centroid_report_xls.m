function climada_EDS_ED_at_centroid_report_xls(EDS,xls_file,sheet,varargin)
% climada
% NAME:
%   climada_EDS_ED_at_centroid_report
% PURPOSE:
%   Write out ED at centroids for one or multiple EDS structures into .csv
%   file
%   previous call: climada_EDS_calc
% CALLING SEQUENCE:
%   climada_EDS_ED_at_centroid_report(EDS,entity,xls_file)
% EXAMPLE:
%   climada_EDS_ED_at_centroid_report(climada_EDS_calc(climada_entity_read), climada_entity_read)
% INPUTS:
%   EDS: either an event damage set, as e.g. returned by climada_EDS_calc or
%       a file containing such a structure
%       SPECIAL: we also accept a structure which contains an EDS, like
%       measures_impact.EDS
%       if EDS has the field annotation_name, the legend will show this
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
%   xls_file: filename (and path) to save the report to (as .xls), if
%       empty, prompted for
% OUTPUTS:
%   none, report file written as .xls
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150430, init
% Gilles Stassen, gillesstassen@hotmail.com, 20150625, generalisation with varargin
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('EDS'     ,'var'),    EDS     =[];	end
if ~exist('xls_file','var'),    xls_file='';    end
if ~exist('sheet'   ,'var'),    sheet   ='';    end

% PARAMETERS
% prompt for EDS if not given
if isempty(EDS) % local GUI
    EDS=[climada_global.data_dir filesep 'results' filesep '*.mat'];
    [filename, pathname] = uigetfile(EDS, 'Select EDS:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        EDS=fullfile(pathname,filename);
    end
end

% load the entity, if a filename has been passed
if ~isstruct(EDS)
    EDS_file=EDS;EDS=[];
    load(EDS_file);
end

if exist('measures_impact','var') % if a results file is loaded
    EDS=measures_impact.EDS;
end

if isfield(EDS,'EDS')
    EDS      = EDS.EDS;
end

% check if field 'ED_at_centroid' exists
if ~isfield(EDS,'ED_at_centroid')
    fprintf('Field ''ED_at_centroid'' not provided in EDS structure.\n')
    return
else
    EDS_to_print = [];
    for e_i = 1:length(EDS)
        if length(EDS(e_i).ED_at_centroid) ~= length(EDS(1).ED_at_centroid)
            cprintf([1 0.5 0],'WARNING: length of data field ED_at_centroid of EDS %i incompatible with EDS(1), skipped\n',e_i)
        else
            EDS_to_print = [EDS_to_print e_i];
        end
    end
end
            
flds = {};
if ~isempty(varargin)
    for var_i = 1:length(varargin)
        if isfield(EDS,varargin{var_i})
            ok = [];
            for e_i = EDS_to_print
                if length(EDS(e_i).(varargin{var_i})) == length(EDS(1).ED_at_centroid)
                    ok = [ok 1];
                else
                    ok = [ok 0];
                    cprintf([1 0.5 0],'WARNING: length of data field %s of EDS %i incompatible with EDS(1), skipped\n',varargin{var_i},e_i)
                end
            end
            if all(ok),  flds{end+1} = varargin{var_i};   end
        else
            cprintf([1 0.5 0],'WARNING: fieldname %s not found\n',varargin{var_i})
        end
    end
end

if isempty(xls_file)
    xls_file=[climada_global.data_dir filesep 'results' filesep 'ED_report.xls'];
    [filename, pathname] = uiputfile(xls_file, 'Save ED at centroid report as:');
    if isequal(filename,0) || isequal(pathname,0)
        xls_file='';
    else
        xls_file=fullfile(pathname,filename);
    end
end


% write data into matrix, which will be outputted to xls
msgstr = '& '; %init
for fld_i = 1:length(flds)
    msgstr = [msgstr strrep(flds{fld_i},'_',' ') ', '];
end
fprintf('writing ED at centroid %s\b\b to xls... ',msgstr)

matr          = cell(length(EDS(1).ED_at_centroid)+4,numel(EDS)*2+2);
matr{4,1}     = 'X';
matr{4,2}     = 'Y';
matr(5:end,1) = num2cell(EDS(1).assets.X);
matr(5:end,2) = num2cell(EDS(1).assets.Y);

static_col = 2;

% special case for Barisal where we have two additional variables to
% describee the centroids
if exist(EDS(1).assets.filename,'file')
    load(EDS(1).assets.filename)
    if isfield(entity(1).assets,'Att_100x100_Cell_Code')
        matr{4,static_col+1}     = '100x100 Cell code';
        matr(5:end,static_col+1) = entity(1).assets.Att_100x100_Cell_Code;
        static_col = static_col +1;
    end
    if isfield(entity(1).assets,'Ward_Nr')
        matr{4,static_col+1}     = 'Ward no';
        matr(5:end,static_col+1) = num2cell(entity(1).assets.Ward_Nr);
        static_col = static_col +1;
    end
    if isfield(entity(1).assets,'Category')
        matr{4,static_col+1}     = 'Category';
        matr(5:end,static_col+1) = entity(1).assets.Category;
        static_col = static_col +1;
    end
end

single_Value_col = 1; sim_ndx = ones(1,length(EDS(1).assets.filename)); % init
for e_i = 1:length(EDS)
    if any(EDS(e_i).assets.Value ~= EDS(1).assets.Value),    single_Value_col = 0;   end
%     sim_ndx = sim_ndx & (EDS(e_i).assets.filename == EDS(1).assets.filename);
end

if single_Value_col    
    n_cols  = 1+length(flds); 
    col_1   = 1;
    
    [~,entity_name] = fileparts(EDS(1).assets.filename(sim_ndx));
    entity_name = strrep(entity_name,'_',' ');
    
    matr{1,     static_col+1} = sprintf('Total value');
    matr{2,     static_col+1} = sum(EDS(e_i).assets.Value);
    matr{4,     static_col+1} = sprintf('Total exposure value %s',entity_name);
    matr(5:end, static_col+1) = num2cell(entity(1).assets.Value);
    static_col = static_col +1;
else
    n_cols  = 2+length(flds);
    col_1   = 2;
end

for e_i = 1:numel(EDS)
    [~,entity_name] = fileparts(EDS(e_i).assets.filename);
    entity_name = strrep(entity_name,'_',' ');
    
    if ~single_Value_col
        matr{1,     (e_i-1)*n_cols+1 +static_col} = sprintf('Total asset value');
        matr{2,     (e_i-1)*n_cols+1 +static_col} = sum(EDS(e_i).assets.Value);
        matr{4,     (e_i-1)*n_cols+1 +static_col} = sprintf('Total asset value %s', entity_name);
        matr(5:end, (e_i-1)*n_cols+1 +static_col) = num2cell(EDS(e_i).assets.Value);
    end
    
    matr{1,     (e_i-1)*n_cols+col_1 +static_col} = sprintf('Total damage (absolute; %%age of TAV)');
    matr{2,     (e_i-1)*n_cols+col_1 +static_col} = sum(EDS(e_i).ED_at_centroid);
    matr{3,     (e_i-1)*n_cols+col_1 +static_col} = sum(EDS(e_i).ED_at_centroid)/sum(EDS(e_i).assets.Value);
    matr{4,     (e_i-1)*n_cols+col_1 +static_col} = sprintf('AED %s',strrep(EDS(e_i).annotation_name,'_',' ')); 
    matr(5:end, (e_i-1)*n_cols+col_1 +static_col) = num2cell(EDS(e_i).ED_at_centroid);

    if isempty(length(flds)),continue; end
    for fld_i =1:length(flds)
        matr{1,     (e_i-1)*n_cols+col_1+fld_i+static_col} = sprintf(['Total ' strrep(flds{fld_i},'_',' ') ' (absolute; %%age of AED)']);
        matr{2,     (e_i-1)*n_cols+col_1+fld_i+static_col} = sum(EDS(e_i).(flds{fld_i}));
        matr{3,     (e_i-1)*n_cols+col_1+fld_i+static_col} = sum(EDS(e_i).(flds{fld_i}))/sum(EDS(end).ED_at_centroid);
        matr{4,     (e_i-1)*n_cols+col_1+fld_i+static_col} = sprintf('%s %s',strrep(flds{fld_i},'_',' '),EDS(e_i).annotation_name);
        matr(5:end, (e_i-1)*n_cols+col_1+fld_i+static_col) = num2cell(EDS(e_i).(flds{fld_i}));
    end
end


warning('off','MATLAB:xlswrite:AddSheet'); % suppress warning message
try
    xlswrite(xls_file,matr,sheet)
catch
    % probably too large for old excel, try writing to .xlsx instead
    try
        xlsx_file = [xls_file 'x'];
        xlswrite(xlsx_file,matr,sheet)
    catch
        % probably too large for new excel, write to textfile instead
        cprintf([1 0 0],'FAILED\n')
        fprintf('attempting to write to text file instead... ')
        txt_file = strrep(xlsx_file,'.xlsx','.txt');
        writetable(cell2table(matr),txt_file)
        fclose all;
    end
end

fprintf('done\n')
if single_Value_col
    cprintf([0 0 1],'NOTE: equivalent asset values for each EDS, truncated into single column\n')
end
fprintf('report written to sheet %s of %s\n',sheet,xls_file);

return
