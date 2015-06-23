
% create Damage-Frequency-Graph for Barisal, showing one DFC per peril

%% Directories
barisal_data_dir= ['\\CHRB1065.CORP.GWPNET.COM\homes\X\S3BXXW\Documents\lea\climada_git\climada_modules\barisal_demo\data'];
% barisal_data_dir= [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
entities_dir    = [barisal_data_dir filesep 'entities'];
hazards_dir     = [barisal_data_dir filesep 'hazards'];
results_dir     = [barisal_data_dir filesep 'results'];

% 
EDS_file= 'EDS_230615';
load([results_dir filesep EDS_file])

%% create DFC graph with different perils
close all
climada_figuresize(0.5,0.9);
[fig,legend_str,return_period,sorted_damage] = climada_EDS_DFC(EDS([3 4 1 2 5]));
print(fig,'-dpdf',[results_dir filesep 'DFC_2014.pdf'])
    
%% create DFC report (xls) with different perils
EDS_no = 15;
matr          = cell(29+2,EDS_no*2);
for i = 1:EDS_no
    matr{1,(i-1)*2+1} = EDS(i).annotation_name;
    matr{2,(i-1)*2+1} = 'Return period (years)';
    matr{2,(i-1)*2+2} = 'Damage (BDT)';
    [~,~,return_period,sorted_damage] = climada_EDS_DFC(EDS(i));
    close all
    no = numel(return_period);
    matr(3:3+no-1,(i-1)*2+1) = num2cell(return_period);
    matr(3:3+no-1,(i-1)*2+2) = num2cell(sorted_damage);
end
xls_file = [results_dir filesep 'Damage_frequency_curve.xls'];
xlswrite(xls_file,matr)


%%

