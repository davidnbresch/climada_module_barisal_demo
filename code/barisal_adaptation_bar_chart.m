

% barisal measures adaptation bar chart
% ------------------------------------
% create adaptation bar chart based on measures information from table 5.3
% from Vulnerability Analysis Report
% ------------------------------------
% Lea Mueller, muellele@gmail.com, 20151125, rename to climada_adaptation_bar_chart from climada_adaptation_bar_chart_v2


% set directory
barisal_data_dir= [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
results_dir     = [barisal_data_dir filesep 'results' filesep '20150916_new_runs'];
return

% read input from excel (table 5.3)
excel_file = [climada_global.project_dir filesep 'results' filesep '20150916_new_runs' filesep 'measures' filesep 'measures_impact.xlsx'];
measures = climada_xlsread(0,excel_file,'measures',0);


% create measures_impact structure and set values according to xls-input

% today (cost capex+opex, benefit 2015)
measures_impact.Value_unit = 'Tk. Crore';
measures_impact.title_str = 'Benefit 2015';
% measures_impact.NPV_total_climate_risk = 100;
measures_impact.NPV_total_climate_risk = 90;
measures_impact.benefit = measures.benefit_2015;
measures_impact.cb_ratio = 1./measures.bc_ratio;
measures_impact.measures.name = measures.name;
measures_impact.measures.cost = measures.cost_capex_opex;

% reference (cost capex benefit 2050)
% measures_impact_reference = measures_impact;
% measures_impact_reference.benefit = measures.benefit_2050;
% measures_impact_reference.measures.cost = measures.cost_capex;
measures_impact(2) = measures_impact;
measures_impact(2).benefit = measures.benefit_2050;
measures_impact(2).measures.cost = measures.cost_capex;
measures_impact(2).title_str = 'Benefit 2050';

save([results_dir filesep 'measures' filesep 'measures_impact_BDT'],'measures_impact');

% create adaptation bar chart
benefit_str = 'Risk reduction in AED (%)';
scale_benefit = 2.5;
tcr_off = 1;
sort_benefit = 1;
fig = climada_adaptation_bar_chart(measures_impact,sort_benefit,scale_benefit,benefit_str,'south',tcr_off);
print(fig,'-dpdf',[results_dir filesep sprintf('Barisal_measures_proposal.pdf')]);
    
% sort_benefit = 0;
% fig = climada_adaptation_bar_chart(measures_impact,measures_impact_reference,sort_benefit,'',sort_benefit,scale_benefit,benefit_str,'east',tcr_off);
% print(fig,'-dpdf',[results_dir filesep sprintf('Barisal_measures_proposal_sorted.pdf')]);

% sort_benefit = 0;
% fig = climada_adaptation_bar_chart(measures_impact,measures_impact_reference,sort_benefit,'',sort_benefit,scale_benefit,benefit_str,'east',tcr_off);
% print(fig,'-dpdf',[results_dir filesep sprintf('Barisal_measures_proposal_sorted.pdf')]);


