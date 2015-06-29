function measures = climada_measures_construct(measures,n_measures)

if ~exist('measures'    ,'var'),    measures = struct([]);  end
if ~exist('n_measures'  ,'var'),    n_measures = 1;         end

if isempty(measures)
    measures.name={};measures.color={};measures.color_RGB = [];
    measures.cost= [];measures.hazard_intensity_impact=[];
    measures.hazard_high_frequency_cutoff=[];measures.hazard_event_set={};
    measures.MDD_impact_a= [];measures.MDD_impact_b= [];measures.PAA_impact_a= [];
    measures.PAA_impact_b= [];measures.damagefunctions_map={}; measures.entity_file={};
    measures.risk_transfer_attachement = [];measures.risk_transfer_cover = [];
    measures.peril_ID={}; measures.hazard_event_set_operator={}; 
end

if n_measures == 0
    return;
end

if n_measures >=1
    for measure_i = 1:n_measures
        measures.name{end+1}    = ['measure_' num2str(measure_i+length(measures.cost))];

        R = rand; G = rand; B = rand; % random colors for a (pleasant) surprise each time :)
        measures.color{end+1}                           = [num2str(R) ' ' num2str(G) ' ' num2str(B)];
        measures.color_RGB(end+1,:)                     = [R; G; B];
        measures.cost(end+1)                            = [1];
        measures.hazard_intensity_impact(end+1)         = [0];
        measures.hazard_high_frequency_cutoff(end+1)    = [0];
        measures.hazard_event_set{end+1}                = 'nil';
        measures.MDD_impact_a(end+1)                    = [1];
        measures.MDD_impact_b(end+1)                    = [0];
        measures.PAA_impact_a(end+1)                    = [1];
        measures.PAA_impact_b(end+1)                    = [0];
        measures.damagefunctions_map{end+1}             = 'nil';
        measures.risk_transfer_attachement(end+1)       = [0];
        measures.risk_transfer_cover(end+1)             = [0];
        measures.entity_file{end+1}                     = 'nil';
        measures.peril_ID{end+1}                        = '';
        measures.hazard_event_set_operator{end+1}       = '';
    end
end

if all(n_measures < 0)
    rm_measures = sort(abs(unique(n_measures)),'descend');
    flds = fieldnames(measures);
    for measure_i = rm_measures
        total_no_measures = length(measures.name);
        for fld_i = 1:length(flds)
            if length(measures.(flds{fld_i})) == total_no_measures && ~ischar(measures.(flds{fld_i}))
                if strcmp(flds{fld_i},'color_RGB')
                    measures.(flds{fld_i})(measure_i,:) = [];
                else
                    measures.(flds{fld_i})(measure_i) = [];
                end
            end
        end
    end
end                  

measures = climada_measures_encode(measures);