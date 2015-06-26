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
if ~exist('varargin','var'),    varargin='';    end

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
            for e_i = EDS_to_print
                if length(EDS.(varargin{var_i})) == length(EDS(1).ED_at_centroid)
                    flds{end+1} = varargin{var_i};
                else
                    cprintf([1 0.5 0],'WARNING: length of data field %s of EDS %i incompatible with EDS(1), skipped\n',varargin{var_i},e_i)
                end
            end
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
    msgstr = [msgstr flds{fld_i} ', '];
end
fprintf('writing ED at centroid %s\b\b to xls... ',msgstr)

matr          = cell(length(EDS(1).ED_at_centroid)+3,numel(EDS)*2+2);
matr{3,1}     = 'X';
matr{3,2}     = 'Y';
matr(4:end,1) = num2cell(EDS(1).assets.X);
matr(4:end,2) = num2cell(EDS(1).assets.Y);

static_row_no = 2;

% special case for Barisal where we have two additional variables to
% describee the centroids
if exist(EDS(1).assets.filename,'file')
    load(EDS(1).assets.filename)
end
if isfield(entity(1).assets,'Att_100x100_Cell_Code')
    matr{3,end+1}     = '100x100 Cell code';
    matr(4:end,end+1) = entity(1).assets.Att_100x100_Cell_Code;
    static_row_no = static_row_no +1;
end
if isfield(entity(1).assets,'Ward_Nr')
    matr{3,end+1}     = 'Ward no';
    matr(4:end,end+1) = num2cell(entity(1).assets.Ward_Nr);
    static_row_no = static_row_no +1;
end
if isfield(entity(1).assets,'Category')
    matr{3,end+1}     = 'Category';
    matr(4:end,end+1) = entity(1).assets.Category;
    static_row_no = static_row_no +1;
end

n_flds = 2+length(flds);

for e_i = 1:numel(EDS)
    [~,entity_name] = fileparts(EDS(e_i).assets.filename);
    entity_name = strrep(entity_name,'_',' ');
    
    matr{3,(e_i-1)*n_flds+1 +static_row_no} = sprintf('Total exposure value %s', entity_name);
    matr{3,(e_i-1)*n_flds+2 +static_row_no} = sprintf('AED %s',EDS(e_i).annotation_name); 
    
    matr(4:end,(e_i-1)*n_flds+1 +static_row_no) = num2cell(EDS(e_i).assets.Value);
    matr(4:end,(e_i-1)*n_flds+2 +static_row_no) = num2cell(EDS(e_i).ED_at_centroid);

    matr{1,(e_i-1)*n_flds+1 +static_row_no} = sprintf('Total Value');
    matr{2,(e_i-1)*n_flds+1 +static_row_no} = sum(EDS(e_i).assets.Value);

    matr{1,(e_i-1)*n_flds+2 +static_row_no} = sprintf('Total damage');
    matr{2,(e_i-1)*n_flds+2 +static_row_no} = sum(EDS(e_i).ED_at_centroid);
    
    if isempty(length(flds)),continue; end
    for fld_i = 3:length(n_flds)
        matr{3,(e_i-1)*n_flds+fld_i +static_row_no} = sprintf('%s %s',strrep(flds{fld_i},'_',' '),EDS(e_i).annotation_name);
        
        matr(4:end,(e_i-1)*n_flds+fld_i +static_row_no) = num2cell(EDS(e_i).(flds{fld_i}));
        
        matr{1,(e_i-1)*n_flds+fld_i +static_row_no} = sprintf(['Total' strrep(flds{fld_i},'_',' ')]);
        matr{2,(e_i-1)*n_flds+fld_i +static_row_no} = sum(EDS(e_i).(flds{fld_i}));
    end
end
try
    xlswrite(xls_file,matr,sheet)
catch
    % probably too large for old excel, try writing to .xlsx instead
    xlsx_file = [xls_file 'x'];
    xlswrite(xlsx_file,matr,sheet)
end

fprintf('done\n')

fprintf('report written to sheet %s of %s\n',sheet,xls_file);

return
