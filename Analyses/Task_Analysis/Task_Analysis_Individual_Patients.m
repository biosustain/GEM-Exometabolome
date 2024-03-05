initCobraToolbox (0)
changeCobraSolver('gurobi') 


model = readCbModel('sampleMetaOut_EndoRecon1NewBounds.mat');
model_points = model.points; 

tasks_table=readtable('Supp_Tables_2.xls','sheet','Supp_Table_1_Recon1_Filtered');
tasks=unique(tasks_table(:,1),'stable'); 
tasks_name=unique(tasks_table(:,4),'stable'); 
list_of_tasks=unique(join(string(tasks_table{:,2:4}),'_'),'stable');
list_of_tasks_split=split(list_of_tasks,'_');

list_of_patient_models = readtable('List_Of_Patient_Models.xls','sheet','sheet1');
list_of_patient_id = readtable('List_Of_Patient_Models.xls','sheet','sheet2'); % The ID contains the patient's number and group
summary_table=zeros(height(tasks),height(list_of_patient_models));

for j=1:height(list_of_patient_models)
    % Load the jth GEM
    ith_model = readCbModel(list_of_patient_models{j,1}{1});
    for i=1:height(tasks)
        % Ith Function
        ith_Index=find(ismember(tasks_table(:,1),tasks(i,1)));  
        ith_substrate=table2cell(tasks_table(ith_Index,5));
        ith_product=table2cell(tasks_table(ith_Index,8));
        % Testing patient
	    [flux, ~ ,~]=testPathway(ith_model,ith_substrate,ith_product);
        summary_table(i,j) = flux;
    end    
end
writematrix(summary_table, 'Patients_Task_Analysis.xls', 'Sheet', 1);
save('Work_Space-Task_Analyisi.mat')

cgo=clustergram(summary_table,'Standardize','Row', 'Colormap', 'redbluecmap');
%set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'TC-P1', 'TA-P2', 'TC-P3', 'TD-P4', 'TD-P5', 'TD-P6', 'TA-P7', 'TA-P8', 'TB-P9', 'TD-P10', 'TC-P11', 'TD-P12', 'TD-P13', 'TC-P14', 'TD-P15', 'TC-P16', 'TD-P17', 'TC-P18', 'TC-P19', 'TB-P20', 'TB-P21', 'TD-P22', 'TA-P23', 'TD-P24', 'TD-P25', 'TD-P26', 'TD-P27', 'TD-P28', 'TB-P29', 'TC-P30', 'TC-P31', 'TD-P32', 'TD-P33', 'TC-P34', 'TD-P35', 'TB-P36', 'TB-P37', 'TC-P38', 'TC-P39', 'TC-P40', 'TC-P41', 'TC-P42', 'TB-P43', 'TD-P44', 'TC-P45', 'TD-P46', 'TD-P47', 'TB-P48', 'TC-P49', 'TD-P50', 'TB-P51', 'TA-P52', 'TA-P53', 'TD-P54', 'TD-P55', 'TB-P56', 'TD-P57', 'TA-P58', 'TC-P59', 'TD-P60', 'TC-P61', 'TC-P62', 'TC-P63', 'TC-P64', 'TD-P65', 'TD-P66', 'TD-P67', 'TD-P68', 'TD-P69', 'TA-P70', 'TB-P71', 'TD-P72', 'TA-P73', 'TD-P74', 'TD-P75', 'TA-P76', 'TC-P77', 'TD-P78', 'TC-P79', 'TC-P80', 'TD-P81', 'TC-P82', 'TD-P83', 'TD-P84', 'TB-P85', 'TC-P86', 'TD-P87', 'TD-P88', 'TD-P89', 'TA-P90', 'TD-P91', 'TD-P92', 'TC-P93', 'TA-P94', 'TC-P95'}, 'RowLabels', {'Oxidative phosphorylation via NADH-coenzyme Q oxidoreductase (COMPLEX I)', 'Oxidative phosphorylation via succinate-coenzyme Q oxidoreductase (COMPLEX II)', 'Krebs cycle - oxidative decarboxylation of pyruvate', 'Krebs cycle - NADH generation', 'ATP regeneration from glucose (normoxic conditions) - glycolysis + krebs cycle', 'ATP generation from glucose (hypoxic conditions) - glycolysis', 'Reactive oxygen species detoxification (H2O2 to H2O)', 'Presence of the thioredoxin system through the thioredoxin reductase activity', 'Inosine monophosphate synthesis (IMP)', 'Cytidine triphosphate synthesis (CTP)', 'Guanosine triphosphate synthesis (GTP)', 'Uridine triphosphate synthesis (UTP)', 'Deoxyadenosine triphosphate synthesis (dATP)', 'Deoxycytidine triphosphate synthesis (dCTP)', 'Deoxyguanosine triphosphate synthesis (dGTP)', 'Deoxyuridine triphosphate synthesis (dUTP)', 'Deoxythymidine triphosphate synthesis (dTTP)', 'AMP salvage from adenine', 'IMP salvage from hypoxanthine', 'GMP salvage from guanine', '3-Phospho-5-adenylyl sulfate synthesis', 'Degradation of adenine to urate', 'Degradation of guanine to urate', 'Degradation of cytosine', 'Degradation of uracil', 'Gluconeogenesis from pyruvate', 'Gluconeogenesis from Lactate', 'Gluconeogenesis from Glycerol', 'Gluconeogenesis from Alanine', 'Gluconeogenesis from Glutamine', 'Ethanol to acetaldehyde', 'Glucose to lactate conversion', 'Malate to pyruvate conversion', 'Synthesis of fructose-6-phosphate from erythrose-4-phosphate (HMP shunt)', 'Synthesis of ribose-5-phosphate', 'Synthesis of lactose', 'Glycogen biosynthesis', 'Glycogen degradation', 'Fructose degradation (to glucose-3-phosphate)', 'Fructose to glucose conversion (via fructose-6-phosphate)', 'UDP-glucose synthesis', 'UDP-galactose synthesis', 'UDP-glucuronate synthesis', 'GDP-L-fucose synthesis', 'Mannose degradation (to fructose-6-phosphate)', 'GDP-mannose synthesis', 'UDP-N-acetyl D-galactosamine synthesis', 'CMP-N-acetylneuraminate synthesis', 'N-Acetylglucosamine synthesis', 'Glucuronate synthesis (via inositol)', 'Glucuronate synthesis (via udp-glucose)', 'Synthesis of acetone', '(R)-3-Hydroxybutanoate synthesis', 'Synthesis of inositol', 'Inositol as input for glucuronate-xylulose pathway', 'Synthesis of phosphatidylinositol from inositol', 'Conversion of phosphatidyl-1D-myo-inositol to 1D-myo-inositol 1-phosphate', 'Starch degradation', 'Link between glyoxylate metabolism and pentose phosphate pathway (Xylulose to glycolate)', 'Synthesis of methylglyoxal', 'Alanine synthesis', 'Alanine degradation', 'Synthesis of alanine from glutamine', 'Arginine synthesis', 'Arginine degradation', 'Synthesis of arginine from glutamine', 'Synthesis of nitric oxide from arginine', 'Synthesis of aspartate from glutamine', 'Synthesis of creatine from arginine', 'Asparagine synthesis', 'Asparagine degradation', 'Aspartate synthesis', 'Aspartate degradation', 'Conversion of aspartate to arginine', 'Conversion of aspartate to beta-alanine', 'Conversion of asparate to asparagine', 'beta-Alanine synthesis', 'beta-Alanine degradation', 'Beta-alanine to 3-oxopropanoate', 'Cysteine synthesis (need serine and methionine)', 'Cysteine degradation', 'Synthesis of cysteine from cystine', 'Synthesis of taurine from cysteine', 'Glutamate synthesis', 'Glutamate degradation', 'Conversion of glutamate to glutamine', 'Conversion of glutamate to proline', 'Conversion of GABA into succinate', 'Glutamine synthesis', 'Glutamine degradation', 'Glutaminolysis (glutamine to lactate)', 'Glutathionate synthesis', 'Glycine synthesis', 'Glycine degradation', 'Conversion of glycine to pyruvate', 'Histidine degradation', 'Conversion of histidine to histamine', 'Homocysteine synthesis (need methionine)', 'Homocysteine degradation', 'Isoleucine degradation', 'Leucine degradation', 'Conversion of leucine to acetyl-coA', 'Lysine degradation', 'Conversion of lysine to L-Saccharopine', 'Conversion of lysine to L-2-Aminoadipate', 'Methionine degradation', 'S-adenosyl-L-methionine synthesis', 'Ornithine degradation', 'Synthesis of ornithine from glutamine', 'Synthesis of spermidine from ornithine', 'Serine synthesis', 'Serine degradation', 'Phenylalanine degradation', 'Phenylalanine to phenylacetaldehyde', 'Phenylalanine to phenylacetyl-L-glutaminate', 'Phenylalanine to tyrosine', 'Proline synthesis', 'Proline degradation', 'Synthesis of proline from glutamine', 'Threonine degradation', 'Tryptophan degradation', 'Synthesis of anthranilate from tryptophan', 'Synthesis of L-kynurenine from tryptophan', 'Synthesis of N-formylanthranilate from tryptophan', 'Synthesis of quinolinate from tryptophan', 'Synthesis of serotonin from tryptophan', 'Tyrosine synthesis (need phenylalanine)', 'Tyrosine degradation', 'Tyrosine to adrenaline', 'Tyrosine to dopamine', 'Tyrosine to acetoacetate and fumarate', 'Valine degradation', 'Valine to succinyl-coA', 'Hydroxymethylglutaryl-CoA synthesis', 'Cholesterol synthesis', 'Acetoacetate synthesis', 'Mevalonate synthesis', 'Farnesyl-pyrophosphate synthesis', 'Glycerol-3-phosphate synthesis', 'Phosphatidyl-choline synthesis', 'Phosphatidyl-ethanolamine synthesis', 'Phosphatidyl-serine synthesis', 'Phosphatidyl-inositol synthesis', 'Cardiolipin synthesis', 'Triacylglycerol synthesis', 'Sphingomyelin synthesis', 'Ceramide synthesis', 'Palmitate synthesis', 'Palmitate degradation', 'Palmitolate synthesis', 'Palmitolate degradation', 'cis-vaccenic acid synthesis', 'cis-vaccenic acid degradation', 'Elaidate synthesis', 'Elaidate degradation', 'Linolenate degradation', 'Linoleate degradation', 'gamma-Linolenate synthesis', 'gamma-Linolenate degradation', 'Arachidonate synthesis', 'Arachidonate degradation', 'Synthesis of malonyl-coa', 'Synthesis of palmitoyl-CoA', 'Taurochenodeoxycholate synthesis', 'Glycochenodeoxycholate synthesis', 'tauro-cholate synthesis', 'glyco-cholate synthesis', 'Synthesis of thromboxane from arachidonate', 'Synthesis of galactosyl glucosyl ceramide (link with ganglioside metabolism)', 'Synthesis of glucocerebroside', 'Synthesis of globoside (link with globoside metabolism)', 'NAD synthesis from nicotinamide', 'FAD synthesis', 'Synthesis of coenzyme A', 'Tetrahydrofolate synthesis', 'Pyridoxal-phosphate synthesis', 'Heme synthesis', 'Phosphatidyl-inositol to glucosaminyl-acylphosphatidylinositol', 'Glucosaminyl-acylphosphatidylinositoll to deacylated-glycophosphatidylinositol (GPI)-anchored protein', 'Degradation of n2m2nmasn', 'Biosynthesis of m4mpdol-U', 'Biosynthesis of g3m8masn', 'Degradation of s2l2fn2m2masn', 'Biosynthesis of core2 (beta-D-Galactosyl-1,3-(N-acetyl-beta-D-glucosaminyl-1,6)-N-acetyl-D-galactosaminyl-R)', 'Biosynthesis of core4 (NAcetyl-beta-D-glucosaminyl-1,6-(Nacetyl-beta-D-glucosaminyl-1,3)-N-acetyl-D-galactosaminyl-R)', 'Biosynthesis of Tn-antigen (Glycoprotein N-acetyl-D-galactosamine)', 'Keratan sulfate biosynthesis from O-glycan (core 2-linked)', 'Keratan sulfate biosynthesis from O-glycan (core 4-linked)', 'Keratan sulfate degradation', 'Keratan sulfate biosynthesis from N-glycan'})%, 'ColumnLabelsRotate', 45)
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'TC-P1', 'TA-P2', 'TC-P3', 'TD-P4', 'TD-P5', 'TD-P6', 'TA-P7', 'TA-P8', 'TB-P9', 'TD-P10', 'TC-P11', 'TD-P12', 'TD-P13', 'TC-P14', 'TD-P15', 'TC-P16', 'TD-P17', 'TC-P18', 'TC-P19', 'TB-P20', 'TB-P21', 'TD-P22', 'TA-P23', 'TD-P24', 'TD-P25', 'TD-P26', 'TD-P27', 'TD-P28', 'TB-P29', 'TC-P30', 'TC-P31', 'TD-P32', 'TD-P33', 'TC-P34', 'TD-P35', 'TB-P36', 'TB-P37', 'TC-P38', 'TC-P39', 'TC-P40', 'TC-P41', 'TC-P42', 'TB-P43', 'TD-P44', 'TC-P45', 'TD-P46', 'TD-P47', 'TB-P48', 'TC-P49', 'TD-P50', 'TB-P51', 'TA-P52', 'TA-P53', 'TD-P54', 'TD-P55', 'TB-P56', 'TD-P57', 'TA-P58', 'TC-P59', 'TD-P60', 'TC-P61', 'TC-P62', 'TC-P63', 'TC-P64', 'TD-P65', 'TD-P66', 'TD-P67', 'TD-P68', 'TD-P69', 'TA-P70', 'TB-P71', 'TD-P72', 'TA-P73', 'TD-P74', 'TD-P75', 'TA-P76', 'TC-P77', 'TD-P78', 'TC-P79', 'TC-P80', 'TD-P81', 'TC-P82', 'TD-P83', 'TD-P84', 'TB-P85', 'TC-P86', 'TD-P87', 'TD-P88', 'TD-P89', 'TA-P90', 'TD-P91', 'TD-P92', 'TC-P93', 'TA-P94', 'TC-P95'})%, 'ColumnLabelsRotate', 45)
print(gcf, 'Task_Analysis_Row_Clustering_2.svg', '-dsvg'); % Save the figure as an SVG file

cgo=clustergram(summary_table,'Standardize','Column', 'Colormap', 'redbluecmap');
%set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'TC-P1', 'TA-P2', 'TC-P3', 'TD-P4', 'TD-P5', 'TD-P6', 'TA-P7', 'TA-P8', 'TB-P9', 'TD-P10', 'TC-P11', 'TD-P12', 'TD-P13', 'TC-P14', 'TD-P15', 'TC-P16', 'TD-P17', 'TC-P18', 'TC-P19', 'TB-P20', 'TB-P21', 'TD-P22', 'TA-P23', 'TD-P24', 'TD-P25', 'TD-P26', 'TD-P27', 'TD-P28', 'TB-P29', 'TC-P30', 'TC-P31', 'TD-P32', 'TD-P33', 'TC-P34', 'TD-P35', 'TB-P36', 'TB-P37', 'TC-P38', 'TC-P39', 'TC-P40', 'TC-P41', 'TC-P42', 'TB-P43', 'TD-P44', 'TC-P45', 'TD-P46', 'TD-P47', 'TB-P48', 'TC-P49', 'TD-P50', 'TB-P51', 'TA-P52', 'TA-P53', 'TD-P54', 'TD-P55', 'TB-P56', 'TD-P57', 'TA-P58', 'TC-P59', 'TD-P60', 'TC-P61', 'TC-P62', 'TC-P63', 'TC-P64', 'TD-P65', 'TD-P66', 'TD-P67', 'TD-P68', 'TD-P69', 'TA-P70', 'TB-P71', 'TD-P72', 'TA-P73', 'TD-P74', 'TD-P75', 'TA-P76', 'TC-P77', 'TD-P78', 'TC-P79', 'TC-P80', 'TD-P81', 'TC-P82', 'TD-P83', 'TD-P84', 'TB-P85', 'TC-P86', 'TD-P87', 'TD-P88', 'TD-P89', 'TA-P90', 'TD-P91', 'TD-P92', 'TC-P93', 'TA-P94', 'TC-P95'}, 'RowLabels', {'Oxidative phosphorylation via NADH-coenzyme Q oxidoreductase (COMPLEX I)', 'Oxidative phosphorylation via succinate-coenzyme Q oxidoreductase (COMPLEX II)', 'Krebs cycle - oxidative decarboxylation of pyruvate', 'Krebs cycle - NADH generation', 'ATP regeneration from glucose (normoxic conditions) - glycolysis + krebs cycle', 'ATP generation from glucose (hypoxic conditions) - glycolysis', 'Reactive oxygen species detoxification (H2O2 to H2O)', 'Presence of the thioredoxin system through the thioredoxin reductase activity', 'Inosine monophosphate synthesis (IMP)', 'Cytidine triphosphate synthesis (CTP)', 'Guanosine triphosphate synthesis (GTP)', 'Uridine triphosphate synthesis (UTP)', 'Deoxyadenosine triphosphate synthesis (dATP)', 'Deoxycytidine triphosphate synthesis (dCTP)', 'Deoxyguanosine triphosphate synthesis (dGTP)', 'Deoxyuridine triphosphate synthesis (dUTP)', 'Deoxythymidine triphosphate synthesis (dTTP)', 'AMP salvage from adenine', 'IMP salvage from hypoxanthine', 'GMP salvage from guanine', '3-Phospho-5-adenylyl sulfate synthesis', 'Degradation of adenine to urate', 'Degradation of guanine to urate', 'Degradation of cytosine', 'Degradation of uracil', 'Gluconeogenesis from pyruvate', 'Gluconeogenesis from Lactate', 'Gluconeogenesis from Glycerol', 'Gluconeogenesis from Alanine', 'Gluconeogenesis from Glutamine', 'Ethanol to acetaldehyde', 'Glucose to lactate conversion', 'Malate to pyruvate conversion', 'Synthesis of fructose-6-phosphate from erythrose-4-phosphate (HMP shunt)', 'Synthesis of ribose-5-phosphate', 'Synthesis of lactose', 'Glycogen biosynthesis', 'Glycogen degradation', 'Fructose degradation (to glucose-3-phosphate)', 'Fructose to glucose conversion (via fructose-6-phosphate)', 'UDP-glucose synthesis', 'UDP-galactose synthesis', 'UDP-glucuronate synthesis', 'GDP-L-fucose synthesis', 'Mannose degradation (to fructose-6-phosphate)', 'GDP-mannose synthesis', 'UDP-N-acetyl D-galactosamine synthesis', 'CMP-N-acetylneuraminate synthesis', 'N-Acetylglucosamine synthesis', 'Glucuronate synthesis (via inositol)', 'Glucuronate synthesis (via udp-glucose)', 'Synthesis of acetone', '(R)-3-Hydroxybutanoate synthesis', 'Synthesis of inositol', 'Inositol as input for glucuronate-xylulose pathway', 'Synthesis of phosphatidylinositol from inositol', 'Conversion of phosphatidyl-1D-myo-inositol to 1D-myo-inositol 1-phosphate', 'Starch degradation', 'Link between glyoxylate metabolism and pentose phosphate pathway (Xylulose to glycolate)', 'Synthesis of methylglyoxal', 'Alanine synthesis', 'Alanine degradation', 'Synthesis of alanine from glutamine', 'Arginine synthesis', 'Arginine degradation', 'Synthesis of arginine from glutamine', 'Synthesis of nitric oxide from arginine', 'Synthesis of aspartate from glutamine', 'Synthesis of creatine from arginine', 'Asparagine synthesis', 'Asparagine degradation', 'Aspartate synthesis', 'Aspartate degradation', 'Conversion of aspartate to arginine', 'Conversion of aspartate to beta-alanine', 'Conversion of asparate to asparagine', 'beta-Alanine synthesis', 'beta-Alanine degradation', 'Beta-alanine to 3-oxopropanoate', 'Cysteine synthesis (need serine and methionine)', 'Cysteine degradation', 'Synthesis of cysteine from cystine', 'Synthesis of taurine from cysteine', 'Glutamate synthesis', 'Glutamate degradation', 'Conversion of glutamate to glutamine', 'Conversion of glutamate to proline', 'Conversion of GABA into succinate', 'Glutamine synthesis', 'Glutamine degradation', 'Glutaminolysis (glutamine to lactate)', 'Glutathionate synthesis', 'Glycine synthesis', 'Glycine degradation', 'Conversion of glycine to pyruvate', 'Histidine degradation', 'Conversion of histidine to histamine', 'Homocysteine synthesis (need methionine)', 'Homocysteine degradation', 'Isoleucine degradation', 'Leucine degradation', 'Conversion of leucine to acetyl-coA', 'Lysine degradation', 'Conversion of lysine to L-Saccharopine', 'Conversion of lysine to L-2-Aminoadipate', 'Methionine degradation', 'S-adenosyl-L-methionine synthesis', 'Ornithine degradation', 'Synthesis of ornithine from glutamine', 'Synthesis of spermidine from ornithine', 'Serine synthesis', 'Serine degradation', 'Phenylalanine degradation', 'Phenylalanine to phenylacetaldehyde', 'Phenylalanine to phenylacetyl-L-glutaminate', 'Phenylalanine to tyrosine', 'Proline synthesis', 'Proline degradation', 'Synthesis of proline from glutamine', 'Threonine degradation', 'Tryptophan degradation', 'Synthesis of anthranilate from tryptophan', 'Synthesis of L-kynurenine from tryptophan', 'Synthesis of N-formylanthranilate from tryptophan', 'Synthesis of quinolinate from tryptophan', 'Synthesis of serotonin from tryptophan', 'Tyrosine synthesis (need phenylalanine)', 'Tyrosine degradation', 'Tyrosine to adrenaline', 'Tyrosine to dopamine', 'Tyrosine to acetoacetate and fumarate', 'Valine degradation', 'Valine to succinyl-coA', 'Hydroxymethylglutaryl-CoA synthesis', 'Cholesterol synthesis', 'Acetoacetate synthesis', 'Mevalonate synthesis', 'Farnesyl-pyrophosphate synthesis', 'Glycerol-3-phosphate synthesis', 'Phosphatidyl-choline synthesis', 'Phosphatidyl-ethanolamine synthesis', 'Phosphatidyl-serine synthesis', 'Phosphatidyl-inositol synthesis', 'Cardiolipin synthesis', 'Triacylglycerol synthesis', 'Sphingomyelin synthesis', 'Ceramide synthesis', 'Palmitate synthesis', 'Palmitate degradation', 'Palmitolate synthesis', 'Palmitolate degradation', 'cis-vaccenic acid synthesis', 'cis-vaccenic acid degradation', 'Elaidate synthesis', 'Elaidate degradation', 'Linolenate degradation', 'Linoleate degradation', 'gamma-Linolenate synthesis', 'gamma-Linolenate degradation', 'Arachidonate synthesis', 'Arachidonate degradation', 'Synthesis of malonyl-coa', 'Synthesis of palmitoyl-CoA', 'Taurochenodeoxycholate synthesis', 'Glycochenodeoxycholate synthesis', 'tauro-cholate synthesis', 'glyco-cholate synthesis', 'Synthesis of thromboxane from arachidonate', 'Synthesis of galactosyl glucosyl ceramide (link with ganglioside metabolism)', 'Synthesis of glucocerebroside', 'Synthesis of globoside (link with globoside metabolism)', 'NAD synthesis from nicotinamide', 'FAD synthesis', 'Synthesis of coenzyme A', 'Tetrahydrofolate synthesis', 'Pyridoxal-phosphate synthesis', 'Heme synthesis', 'Phosphatidyl-inositol to glucosaminyl-acylphosphatidylinositol', 'Glucosaminyl-acylphosphatidylinositoll to deacylated-glycophosphatidylinositol (GPI)-anchored protein', 'Degradation of n2m2nmasn', 'Biosynthesis of m4mpdol-U', 'Biosynthesis of g3m8masn', 'Degradation of s2l2fn2m2masn', 'Biosynthesis of core2 (beta-D-Galactosyl-1,3-(N-acetyl-beta-D-glucosaminyl-1,6)-N-acetyl-D-galactosaminyl-R)', 'Biosynthesis of core4 (NAcetyl-beta-D-glucosaminyl-1,6-(Nacetyl-beta-D-glucosaminyl-1,3)-N-acetyl-D-galactosaminyl-R)', 'Biosynthesis of Tn-antigen (Glycoprotein N-acetyl-D-galactosamine)', 'Keratan sulfate biosynthesis from O-glycan (core 2-linked)', 'Keratan sulfate biosynthesis from O-glycan (core 4-linked)', 'Keratan sulfate degradation', 'Keratan sulfate biosynthesis from N-glycan'})%, 'ColumnLabelsRotate', 45)
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'TC-P1', 'TA-P2', 'TC-P3', 'TD-P4', 'TD-P5', 'TD-P6', 'TA-P7', 'TA-P8', 'TB-P9', 'TD-P10', 'TC-P11', 'TD-P12', 'TD-P13', 'TC-P14', 'TD-P15', 'TC-P16', 'TD-P17', 'TC-P18', 'TC-P19', 'TB-P20', 'TB-P21', 'TD-P22', 'TA-P23', 'TD-P24', 'TD-P25', 'TD-P26', 'TD-P27', 'TD-P28', 'TB-P29', 'TC-P30', 'TC-P31', 'TD-P32', 'TD-P33', 'TC-P34', 'TD-P35', 'TB-P36', 'TB-P37', 'TC-P38', 'TC-P39', 'TC-P40', 'TC-P41', 'TC-P42', 'TB-P43', 'TD-P44', 'TC-P45', 'TD-P46', 'TD-P47', 'TB-P48', 'TC-P49', 'TD-P50', 'TB-P51', 'TA-P52', 'TA-P53', 'TD-P54', 'TD-P55', 'TB-P56', 'TD-P57', 'TA-P58', 'TC-P59', 'TD-P60', 'TC-P61', 'TC-P62', 'TC-P63', 'TC-P64', 'TD-P65', 'TD-P66', 'TD-P67', 'TD-P68', 'TD-P69', 'TA-P70', 'TB-P71', 'TD-P72', 'TA-P73', 'TD-P74', 'TD-P75', 'TA-P76', 'TC-P77', 'TD-P78', 'TC-P79', 'TC-P80', 'TD-P81', 'TC-P82', 'TD-P83', 'TD-P84', 'TB-P85', 'TC-P86', 'TD-P87', 'TD-P88', 'TD-P89', 'TA-P90', 'TD-P91', 'TD-P92', 'TC-P93', 'TA-P94', 'TC-P95'})%, 'ColumnLabelsRotate', 45)
print(gcf, 'Task_Analysis_Col_Clustering_2.svg', '-dsvg'); % Save the figure as an SVG file

% Indexes for Metabogroups
index_ta = find(strcmp(list_of_patient_id{:,2}, 'TA'));
index_tb = find(strcmp(list_of_patient_id{:,2}, 'TB'));
index_tc = find(strcmp(list_of_patient_id{:,2}, 'TC'));
index_td = find(strcmp(list_of_patient_id{:,2}, 'TD'));

%%% Clustering Analysis
%% Plotting
% TA
summary_table_ta = summary_table(:,index_ta);
column_label_ta = list_of_patient_id{index_ta,1};

cgo=clustergram(summary_table_ta,'Standardize','Row', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_ta)%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Row_Clustering.svg', '-dsvg'); % Save the figure as an SVG file

cgo=clustergram(summary_table_ta,'Standardize','Column', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_ta)%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Column_Clustering.svg', '-dsvg'); % Save the figure as an SVG file

%TB
summary_table_tb = summary_table(:,index_tb);
column_label_tb = list_of_patient_id{index_tb,1};

cgo=clustergram(summary_table_tb,'Standardize','Row', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_tb)%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Row_Clustering.svg', '-dsvg'); % Save the figure as an SVG file

cgo=clustergram(summary_table_tb,'Standardize','Column', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_tb)%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Column_Clustering.svg', '-dsvg'); % Save the figure as an SVG file

%TC
summary_table_tc = summary_table(:,index_tc);
column_label_tc = list_of_patient_id{index_tc,1};

cgo=clustergram(summary_table_tc,'Standardize','Row', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_tc)%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Row_Clustering.svg', '-dsvg'); % Save the figure as an SVG file

cgo=clustergram(summary_table_tc,'Standardize','Column', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_tc)%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Column_Clustering.svg', '-dsvg'); % Save the figure as an SVG file

%TD
summary_table_td = summary_table(:,index_td);
column_label_td = list_of_patient_id{index_td,1};

cgo=clustergram(summary_table_td,'Standardize','Row', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_td)%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Row_Clustering.svg', '-dsvg'); % Save the figure as an SVG file

cgo=clustergram(summary_table_td,'Standardize','Column', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_td)%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Column_Clustering.svg', '-dsvg'); % Save the figure as an SVG file


%%Analysis
task_rank = zeros(height(tasks),7); % 7 = number of Metabogroups +3 (progress ,diff between 1-4 -restrictive and non-restrictive-)

%A
cgo=clustergram(summary_table_ta,'Standardize','Column', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_ta, 'RowLabels', table2array(tasks_name))%, 'ColumnLabelsRotate', 45)
task_rank_a = cgo.RowLabels;

%B
cgo=clustergram(summary_table_tb,'Standardize','Column', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_tb, 'RowLabels', table2array(tasks_name))%, 'ColumnLabelsRotate', 45)
task_rank_b = cgo.RowLabels;

%C
cgo=clustergram(summary_table_tc,'Standardize','Column', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_tc, 'RowLabels', table2array(tasks_name))%, 'ColumnLabelsRotate', 45)
task_rank_c = cgo.RowLabels;

%D
cgo=clustergram(summary_table_td,'Standardize','Column', 'Colormap', 'redbluecmap');
set(cgo,'Linkage','complete','Dendrogram',3, 'ColumnLabels', column_label_td, 'RowLabels', table2array(tasks_name))%, 'ColumnLabelsRotate', 45)
task_rank_d = cgo.RowLabels;

task_size = size(tasks_name);
task_significative_change_restrictive = [];
task_significative_change_no_restrictive = [];


for j=1:task_size(1)
    jth_task = tasks_name{j,:}{1};
    task_rank(j,1) = find(strcmp(jth_task, task_rank_a));
    task_rank(j,2) = find(strcmp(jth_task, task_rank_b));
    task_rank(j,3) = find(strcmp(jth_task, task_rank_c));
    task_rank(j,4) = find(strcmp(jth_task, task_rank_d));
    rank_variation = diff([task_rank(j,1),task_rank(j,2),task_rank(j,3),task_rank(j,4)]);
    % restrictive evaluation
    if all(rank_variation < 0)
        task_rank(j,5) = -1;
    elseif all(rank_variation > 0)
        task_rank(j,5) = 1;
    end
    % less restrictive evaluation
    if all(rank_variation <= 0)
        task_rank(j,6) = -1;
    elseif all(rank_variation >= 0)
        task_rank(j,6) = 1;        
    end
    task_rank(j,7) = task_rank(j,4)-task_rank(j,1);
    if abs(task_rank(j,7)) > 10
        if task_rank(j,5) ~= 0
            task_significative_change_restrictive = [task_significative_change_restrictive, j];
        end
        if task_rank(j,6) ~= 0
            task_significative_change_no_restrictive = [task_significative_change_no_restrictive, j];
        end        
    end
end

% save results
tasks_name_column = tasks_name{:,1}; % Extract the column from tasks_name
tasks_name_table = table(tasks_name_column); % Convert tasks_name_column to a table
task_rank_table = array2table(task_rank); % Convert task_rank to a table
combined_table = [tasks_name_table, task_rank_table];
writetable(combined_table, 'Task_analysis_Groups_2.xlsx');



%%% Clustering analysis only with the relevant tasks

%% Restrictive
cgo_restrictive=clustergram(summary_table(task_significative_change_restrictive,:),'Standardize','Row', 'Colormap', 'redbluecmap');
%set(cgo_restrictive,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'TC-P1', 'TA-P2', 'TC-P3', 'TD-P4', 'TD-P5', 'TD-P6', 'TA-P7', 'TA-P8', 'TB-P9', 'TD-P10', 'TC-P11', 'TD-P12', 'TD-P13', 'TC-P14', 'TD-P15', 'TC-P16', 'TD-P17', 'TC-P18', 'TC-P19', 'TB-P20', 'TB-P21', 'TD-P22', 'TA-P23', 'TD-P24', 'TD-P25', 'TD-P26', 'TD-P27', 'TD-P28', 'TB-P29', 'TC-P30', 'TC-P31', 'TD-P32', 'TD-P33', 'TC-P34', 'TD-P35', 'TB-P36', 'TB-P37', 'TC-P38', 'TC-P39', 'TC-P40', 'TC-P41', 'TC-P42', 'TB-P43', 'TD-P44', 'TC-P45', 'TD-P46', 'TD-P47', 'TB-P48', 'TC-P49', 'TD-P50', 'TB-P51', 'TA-P52', 'TA-P53', 'TD-P54', 'TD-P55', 'TB-P56', 'TD-P57', 'TA-P58', 'TC-P59', 'TD-P60', 'TC-P61', 'TC-P62', 'TC-P63', 'TC-P64', 'TD-P65', 'TD-P66', 'TD-P67', 'TD-P68', 'TD-P69', 'TA-P70', 'TB-P71', 'TD-P72', 'TA-P73', 'TD-P74', 'TD-P75', 'TA-P76', 'TC-P77', 'TD-P78', 'TC-P79', 'TC-P80', 'TD-P81', 'TC-P82', 'TD-P83', 'TD-P84', 'TB-P85', 'TC-P86', 'TD-P87', 'TD-P88', 'TD-P89', 'TA-P90', 'TD-P91', 'TD-P92', 'TC-P93', 'TA-P94', 'TC-P95'}, 'RowLabels', {'Oxidative phosphorylation via NADH-coenzyme Q oxidoreductase (COMPLEX I)', 'Oxidative phosphorylation via succinate-coenzyme Q oxidoreductase (COMPLEX II)', 'Krebs cycle - oxidative decarboxylation of pyruvate', 'Krebs cycle - NADH generation', 'ATP regeneration from glucose (normoxic conditions) - glycolysis + krebs cycle', 'ATP generation from glucose (hypoxic conditions) - glycolysis', 'Reactive oxygen species detoxification (H2O2 to H2O)', 'Presence of the thioredoxin system through the thioredoxin reductase activity', 'Inosine monophosphate synthesis (IMP)', 'Cytidine triphosphate synthesis (CTP)', 'Guanosine triphosphate synthesis (GTP)', 'Uridine triphosphate synthesis (UTP)', 'Deoxyadenosine triphosphate synthesis (dATP)', 'Deoxycytidine triphosphate synthesis (dCTP)', 'Deoxyguanosine triphosphate synthesis (dGTP)', 'Deoxyuridine triphosphate synthesis (dUTP)', 'Deoxythymidine triphosphate synthesis (dTTP)', 'AMP salvage from adenine', 'IMP salvage from hypoxanthine', 'GMP salvage from guanine', '3-Phospho-5-adenylyl sulfate synthesis', 'Degradation of adenine to urate', 'Degradation of guanine to urate', 'Degradation of cytosine', 'Degradation of uracil', 'Gluconeogenesis from pyruvate', 'Gluconeogenesis from Lactate', 'Gluconeogenesis from Glycerol', 'Gluconeogenesis from Alanine', 'Gluconeogenesis from Glutamine', 'Ethanol to acetaldehyde', 'Glucose to lactate conversion', 'Malate to pyruvate conversion', 'Synthesis of fructose-6-phosphate from erythrose-4-phosphate (HMP shunt)', 'Synthesis of ribose-5-phosphate', 'Synthesis of lactose', 'Glycogen biosynthesis', 'Glycogen degradation', 'Fructose degradation (to glucose-3-phosphate)', 'Fructose to glucose conversion (via fructose-6-phosphate)', 'UDP-glucose synthesis', 'UDP-galactose synthesis', 'UDP-glucuronate synthesis', 'GDP-L-fucose synthesis', 'Mannose degradation (to fructose-6-phosphate)', 'GDP-mannose synthesis', 'UDP-N-acetyl D-galactosamine synthesis', 'CMP-N-acetylneuraminate synthesis', 'N-Acetylglucosamine synthesis', 'Glucuronate synthesis (via inositol)', 'Glucuronate synthesis (via udp-glucose)', 'Synthesis of acetone', '(R)-3-Hydroxybutanoate synthesis', 'Synthesis of inositol', 'Inositol as input for glucuronate-xylulose pathway', 'Synthesis of phosphatidylinositol from inositol', 'Conversion of phosphatidyl-1D-myo-inositol to 1D-myo-inositol 1-phosphate', 'Starch degradation', 'Link between glyoxylate metabolism and pentose phosphate pathway (Xylulose to glycolate)', 'Synthesis of methylglyoxal', 'Alanine synthesis', 'Alanine degradation', 'Synthesis of alanine from glutamine', 'Arginine synthesis', 'Arginine degradation', 'Synthesis of arginine from glutamine', 'Synthesis of nitric oxide from arginine', 'Synthesis of aspartate from glutamine', 'Synthesis of creatine from arginine', 'Asparagine synthesis', 'Asparagine degradation', 'Aspartate synthesis', 'Aspartate degradation', 'Conversion of aspartate to arginine', 'Conversion of aspartate to beta-alanine', 'Conversion of asparate to asparagine', 'beta-Alanine synthesis', 'beta-Alanine degradation', 'Beta-alanine to 3-oxopropanoate', 'Cysteine synthesis (need serine and methionine)', 'Cysteine degradation', 'Synthesis of cysteine from cystine', 'Synthesis of taurine from cysteine', 'Glutamate synthesis', 'Glutamate degradation', 'Conversion of glutamate to glutamine', 'Conversion of glutamate to proline', 'Conversion of GABA into succinate', 'Glutamine synthesis', 'Glutamine degradation', 'Glutaminolysis (glutamine to lactate)', 'Glutathionate synthesis', 'Glycine synthesis', 'Glycine degradation', 'Conversion of glycine to pyruvate', 'Histidine degradation', 'Conversion of histidine to histamine', 'Homocysteine synthesis (need methionine)', 'Homocysteine degradation', 'Isoleucine degradation', 'Leucine degradation', 'Conversion of leucine to acetyl-coA', 'Lysine degradation', 'Conversion of lysine to L-Saccharopine', 'Conversion of lysine to L-2-Aminoadipate', 'Methionine degradation', 'S-adenosyl-L-methionine synthesis', 'Ornithine degradation', 'Synthesis of ornithine from glutamine', 'Synthesis of spermidine from ornithine', 'Serine synthesis', 'Serine degradation', 'Phenylalanine degradation', 'Phenylalanine to phenylacetaldehyde', 'Phenylalanine to phenylacetyl-L-glutaminate', 'Phenylalanine to tyrosine', 'Proline synthesis', 'Proline degradation', 'Synthesis of proline from glutamine', 'Threonine degradation', 'Tryptophan degradation', 'Synthesis of anthranilate from tryptophan', 'Synthesis of L-kynurenine from tryptophan', 'Synthesis of N-formylanthranilate from tryptophan', 'Synthesis of quinolinate from tryptophan', 'Synthesis of serotonin from tryptophan', 'Tyrosine synthesis (need phenylalanine)', 'Tyrosine degradation', 'Tyrosine to adrenaline', 'Tyrosine to dopamine', 'Tyrosine to acetoacetate and fumarate', 'Valine degradation', 'Valine to succinyl-coA', 'Hydroxymethylglutaryl-CoA synthesis', 'Cholesterol synthesis', 'Acetoacetate synthesis', 'Mevalonate synthesis', 'Farnesyl-pyrophosphate synthesis', 'Glycerol-3-phosphate synthesis', 'Phosphatidyl-choline synthesis', 'Phosphatidyl-ethanolamine synthesis', 'Phosphatidyl-serine synthesis', 'Phosphatidyl-inositol synthesis', 'Cardiolipin synthesis', 'Triacylglycerol synthesis', 'Sphingomyelin synthesis', 'Ceramide synthesis', 'Palmitate synthesis', 'Palmitate degradation', 'Palmitolate synthesis', 'Palmitolate degradation', 'cis-vaccenic acid synthesis', 'cis-vaccenic acid degradation', 'Elaidate synthesis', 'Elaidate degradation', 'Linolenate degradation', 'Linoleate degradation', 'gamma-Linolenate synthesis', 'gamma-Linolenate degradation', 'Arachidonate synthesis', 'Arachidonate degradation', 'Synthesis of malonyl-coa', 'Synthesis of palmitoyl-CoA', 'Taurochenodeoxycholate synthesis', 'Glycochenodeoxycholate synthesis', 'tauro-cholate synthesis', 'glyco-cholate synthesis', 'Synthesis of thromboxane from arachidonate', 'Synthesis of galactosyl glucosyl ceramide (link with ganglioside metabolism)', 'Synthesis of glucocerebroside', 'Synthesis of globoside (link with globoside metabolism)', 'NAD synthesis from nicotinamide', 'FAD synthesis', 'Synthesis of coenzyme A', 'Tetrahydrofolate synthesis', 'Pyridoxal-phosphate synthesis', 'Heme synthesis', 'Phosphatidyl-inositol to glucosaminyl-acylphosphatidylinositol', 'Glucosaminyl-acylphosphatidylinositoll to deacylated-glycophosphatidylinositol (GPI)-anchored protein', 'Degradation of n2m2nmasn', 'Biosynthesis of m4mpdol-U', 'Biosynthesis of g3m8masn', 'Degradation of s2l2fn2m2masn', 'Biosynthesis of core2 (beta-D-Galactosyl-1,3-(N-acetyl-beta-D-glucosaminyl-1,6)-N-acetyl-D-galactosaminyl-R)', 'Biosynthesis of core4 (NAcetyl-beta-D-glucosaminyl-1,6-(Nacetyl-beta-D-glucosaminyl-1,3)-N-acetyl-D-galactosaminyl-R)', 'Biosynthesis of Tn-antigen (Glycoprotein N-acetyl-D-galactosamine)', 'Keratan sulfate biosynthesis from O-glycan (core 2-linked)', 'Keratan sulfate biosynthesis from O-glycan (core 4-linked)', 'Keratan sulfate degradation', 'Keratan sulfate biosynthesis from N-glycan'})%, 'ColumnLabelsRotate', 45)
set(cgo_restrictive,'Linkage','complete','Dendrogram',5, 'ColumnLabels', {'TC-P1', 'TA-P2', 'TC-P3', 'TD-P4', 'TD-P5', 'TD-P6', 'TA-P7', 'TA-P8', 'TB-P9', 'TD-P10', 'TC-P11', 'TD-P12', 'TD-P13', 'TC-P14', 'TD-P15', 'TC-P16', 'TD-P17', 'TC-P18', 'TC-P19', 'TB-P20', 'TB-P21', 'TD-P22', 'TA-P23', 'TD-P24', 'TD-P25', 'TD-P26', 'TD-P27', 'TD-P28', 'TB-P29', 'TC-P30', 'TC-P31', 'TD-P32', 'TD-P33', 'TC-P34', 'TD-P35', 'TB-P36', 'TB-P37', 'TC-P38', 'TC-P39', 'TC-P40', 'TC-P41', 'TC-P42', 'TB-P43', 'TD-P44', 'TC-P45', 'TD-P46', 'TD-P47', 'TB-P48', 'TC-P49', 'TD-P50', 'TB-P51', 'TA-P52', 'TA-P53', 'TD-P54', 'TD-P55', 'TB-P56', 'TD-P57', 'TA-P58', 'TC-P59', 'TD-P60', 'TC-P61', 'TC-P62', 'TC-P63', 'TC-P64', 'TD-P65', 'TD-P66', 'TD-P67', 'TD-P68', 'TD-P69', 'TA-P70', 'TB-P71', 'TD-P72', 'TA-P73', 'TD-P74', 'TD-P75', 'TA-P76', 'TC-P77', 'TD-P78', 'TC-P79', 'TC-P80', 'TD-P81', 'TC-P82', 'TD-P83', 'TD-P84', 'TB-P85', 'TC-P86', 'TD-P87', 'TD-P88', 'TD-P89', 'TA-P90', 'TD-P91', 'TD-P92', 'TC-P93', 'TA-P94', 'TC-P95'}, 'RowLabels', table2cell(tasks_name_table(task_significative_change_restrictive,:)))%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Row_Clustering_2.svg', '-dsvg'); % Save the figure as an SVG file

%% Non-restrictive
cgo_no_restrictive=clustergram(summary_table(task_significative_change_no_restrictive,:),'Standardize','Row', 'Colormap', 'redbluecmap');
%set(cgo_restrictive,'Linkage','complete','Dendrogram',3, 'ColumnLabels', {'TC-P1', 'TA-P2', 'TC-P3', 'TD-P4', 'TD-P5', 'TD-P6', 'TA-P7', 'TA-P8', 'TB-P9', 'TD-P10', 'TC-P11', 'TD-P12', 'TD-P13', 'TC-P14', 'TD-P15', 'TC-P16', 'TD-P17', 'TC-P18', 'TC-P19', 'TB-P20', 'TB-P21', 'TD-P22', 'TA-P23', 'TD-P24', 'TD-P25', 'TD-P26', 'TD-P27', 'TD-P28', 'TB-P29', 'TC-P30', 'TC-P31', 'TD-P32', 'TD-P33', 'TC-P34', 'TD-P35', 'TB-P36', 'TB-P37', 'TC-P38', 'TC-P39', 'TC-P40', 'TC-P41', 'TC-P42', 'TB-P43', 'TD-P44', 'TC-P45', 'TD-P46', 'TD-P47', 'TB-P48', 'TC-P49', 'TD-P50', 'TB-P51', 'TA-P52', 'TA-P53', 'TD-P54', 'TD-P55', 'TB-P56', 'TD-P57', 'TA-P58', 'TC-P59', 'TD-P60', 'TC-P61', 'TC-P62', 'TC-P63', 'TC-P64', 'TD-P65', 'TD-P66', 'TD-P67', 'TD-P68', 'TD-P69', 'TA-P70', 'TB-P71', 'TD-P72', 'TA-P73', 'TD-P74', 'TD-P75', 'TA-P76', 'TC-P77', 'TD-P78', 'TC-P79', 'TC-P80', 'TD-P81', 'TC-P82', 'TD-P83', 'TD-P84', 'TB-P85', 'TC-P86', 'TD-P87', 'TD-P88', 'TD-P89', 'TA-P90', 'TD-P91', 'TD-P92', 'TC-P93', 'TA-P94', 'TC-P95'}, 'RowLabels', {'Oxidative phosphorylation via NADH-coenzyme Q oxidoreductase (COMPLEX I)', 'Oxidative phosphorylation via succinate-coenzyme Q oxidoreductase (COMPLEX II)', 'Krebs cycle - oxidative decarboxylation of pyruvate', 'Krebs cycle - NADH generation', 'ATP regeneration from glucose (normoxic conditions) - glycolysis + krebs cycle', 'ATP generation from glucose (hypoxic conditions) - glycolysis', 'Reactive oxygen species detoxification (H2O2 to H2O)', 'Presence of the thioredoxin system through the thioredoxin reductase activity', 'Inosine monophosphate synthesis (IMP)', 'Cytidine triphosphate synthesis (CTP)', 'Guanosine triphosphate synthesis (GTP)', 'Uridine triphosphate synthesis (UTP)', 'Deoxyadenosine triphosphate synthesis (dATP)', 'Deoxycytidine triphosphate synthesis (dCTP)', 'Deoxyguanosine triphosphate synthesis (dGTP)', 'Deoxyuridine triphosphate synthesis (dUTP)', 'Deoxythymidine triphosphate synthesis (dTTP)', 'AMP salvage from adenine', 'IMP salvage from hypoxanthine', 'GMP salvage from guanine', '3-Phospho-5-adenylyl sulfate synthesis', 'Degradation of adenine to urate', 'Degradation of guanine to urate', 'Degradation of cytosine', 'Degradation of uracil', 'Gluconeogenesis from pyruvate', 'Gluconeogenesis from Lactate', 'Gluconeogenesis from Glycerol', 'Gluconeogenesis from Alanine', 'Gluconeogenesis from Glutamine', 'Ethanol to acetaldehyde', 'Glucose to lactate conversion', 'Malate to pyruvate conversion', 'Synthesis of fructose-6-phosphate from erythrose-4-phosphate (HMP shunt)', 'Synthesis of ribose-5-phosphate', 'Synthesis of lactose', 'Glycogen biosynthesis', 'Glycogen degradation', 'Fructose degradation (to glucose-3-phosphate)', 'Fructose to glucose conversion (via fructose-6-phosphate)', 'UDP-glucose synthesis', 'UDP-galactose synthesis', 'UDP-glucuronate synthesis', 'GDP-L-fucose synthesis', 'Mannose degradation (to fructose-6-phosphate)', 'GDP-mannose synthesis', 'UDP-N-acetyl D-galactosamine synthesis', 'CMP-N-acetylneuraminate synthesis', 'N-Acetylglucosamine synthesis', 'Glucuronate synthesis (via inositol)', 'Glucuronate synthesis (via udp-glucose)', 'Synthesis of acetone', '(R)-3-Hydroxybutanoate synthesis', 'Synthesis of inositol', 'Inositol as input for glucuronate-xylulose pathway', 'Synthesis of phosphatidylinositol from inositol', 'Conversion of phosphatidyl-1D-myo-inositol to 1D-myo-inositol 1-phosphate', 'Starch degradation', 'Link between glyoxylate metabolism and pentose phosphate pathway (Xylulose to glycolate)', 'Synthesis of methylglyoxal', 'Alanine synthesis', 'Alanine degradation', 'Synthesis of alanine from glutamine', 'Arginine synthesis', 'Arginine degradation', 'Synthesis of arginine from glutamine', 'Synthesis of nitric oxide from arginine', 'Synthesis of aspartate from glutamine', 'Synthesis of creatine from arginine', 'Asparagine synthesis', 'Asparagine degradation', 'Aspartate synthesis', 'Aspartate degradation', 'Conversion of aspartate to arginine', 'Conversion of aspartate to beta-alanine', 'Conversion of asparate to asparagine', 'beta-Alanine synthesis', 'beta-Alanine degradation', 'Beta-alanine to 3-oxopropanoate', 'Cysteine synthesis (need serine and methionine)', 'Cysteine degradation', 'Synthesis of cysteine from cystine', 'Synthesis of taurine from cysteine', 'Glutamate synthesis', 'Glutamate degradation', 'Conversion of glutamate to glutamine', 'Conversion of glutamate to proline', 'Conversion of GABA into succinate', 'Glutamine synthesis', 'Glutamine degradation', 'Glutaminolysis (glutamine to lactate)', 'Glutathionate synthesis', 'Glycine synthesis', 'Glycine degradation', 'Conversion of glycine to pyruvate', 'Histidine degradation', 'Conversion of histidine to histamine', 'Homocysteine synthesis (need methionine)', 'Homocysteine degradation', 'Isoleucine degradation', 'Leucine degradation', 'Conversion of leucine to acetyl-coA', 'Lysine degradation', 'Conversion of lysine to L-Saccharopine', 'Conversion of lysine to L-2-Aminoadipate', 'Methionine degradation', 'S-adenosyl-L-methionine synthesis', 'Ornithine degradation', 'Synthesis of ornithine from glutamine', 'Synthesis of spermidine from ornithine', 'Serine synthesis', 'Serine degradation', 'Phenylalanine degradation', 'Phenylalanine to phenylacetaldehyde', 'Phenylalanine to phenylacetyl-L-glutaminate', 'Phenylalanine to tyrosine', 'Proline synthesis', 'Proline degradation', 'Synthesis of proline from glutamine', 'Threonine degradation', 'Tryptophan degradation', 'Synthesis of anthranilate from tryptophan', 'Synthesis of L-kynurenine from tryptophan', 'Synthesis of N-formylanthranilate from tryptophan', 'Synthesis of quinolinate from tryptophan', 'Synthesis of serotonin from tryptophan', 'Tyrosine synthesis (need phenylalanine)', 'Tyrosine degradation', 'Tyrosine to adrenaline', 'Tyrosine to dopamine', 'Tyrosine to acetoacetate and fumarate', 'Valine degradation', 'Valine to succinyl-coA', 'Hydroxymethylglutaryl-CoA synthesis', 'Cholesterol synthesis', 'Acetoacetate synthesis', 'Mevalonate synthesis', 'Farnesyl-pyrophosphate synthesis', 'Glycerol-3-phosphate synthesis', 'Phosphatidyl-choline synthesis', 'Phosphatidyl-ethanolamine synthesis', 'Phosphatidyl-serine synthesis', 'Phosphatidyl-inositol synthesis', 'Cardiolipin synthesis', 'Triacylglycerol synthesis', 'Sphingomyelin synthesis', 'Ceramide synthesis', 'Palmitate synthesis', 'Palmitate degradation', 'Palmitolate synthesis', 'Palmitolate degradation', 'cis-vaccenic acid synthesis', 'cis-vaccenic acid degradation', 'Elaidate synthesis', 'Elaidate degradation', 'Linolenate degradation', 'Linoleate degradation', 'gamma-Linolenate synthesis', 'gamma-Linolenate degradation', 'Arachidonate synthesis', 'Arachidonate degradation', 'Synthesis of malonyl-coa', 'Synthesis of palmitoyl-CoA', 'Taurochenodeoxycholate synthesis', 'Glycochenodeoxycholate synthesis', 'tauro-cholate synthesis', 'glyco-cholate synthesis', 'Synthesis of thromboxane from arachidonate', 'Synthesis of galactosyl glucosyl ceramide (link with ganglioside metabolism)', 'Synthesis of glucocerebroside', 'Synthesis of globoside (link with globoside metabolism)', 'NAD synthesis from nicotinamide', 'FAD synthesis', 'Synthesis of coenzyme A', 'Tetrahydrofolate synthesis', 'Pyridoxal-phosphate synthesis', 'Heme synthesis', 'Phosphatidyl-inositol to glucosaminyl-acylphosphatidylinositol', 'Glucosaminyl-acylphosphatidylinositoll to deacylated-glycophosphatidylinositol (GPI)-anchored protein', 'Degradation of n2m2nmasn', 'Biosynthesis of m4mpdol-U', 'Biosynthesis of g3m8masn', 'Degradation of s2l2fn2m2masn', 'Biosynthesis of core2 (beta-D-Galactosyl-1,3-(N-acetyl-beta-D-glucosaminyl-1,6)-N-acetyl-D-galactosaminyl-R)', 'Biosynthesis of core4 (NAcetyl-beta-D-glucosaminyl-1,6-(Nacetyl-beta-D-glucosaminyl-1,3)-N-acetyl-D-galactosaminyl-R)', 'Biosynthesis of Tn-antigen (Glycoprotein N-acetyl-D-galactosamine)', 'Keratan sulfate biosynthesis from O-glycan (core 2-linked)', 'Keratan sulfate biosynthesis from O-glycan (core 4-linked)', 'Keratan sulfate degradation', 'Keratan sulfate biosynthesis from N-glycan'})%, 'ColumnLabelsRotate', 45)
set(cgo_no_restrictive,'Linkage','complete','Dendrogram',7, 'ColumnLabels', {'TC-P1', 'TA-P2', 'TC-P3', 'TD-P4', 'TD-P5', 'TD-P6', 'TA-P7', 'TA-P8', 'TB-P9', 'TD-P10', 'TC-P11', 'TD-P12', 'TD-P13', 'TC-P14', 'TD-P15', 'TC-P16', 'TD-P17', 'TC-P18', 'TC-P19', 'TB-P20', 'TB-P21', 'TD-P22', 'TA-P23', 'TD-P24', 'TD-P25', 'TD-P26', 'TD-P27', 'TD-P28', 'TB-P29', 'TC-P30', 'TC-P31', 'TD-P32', 'TD-P33', 'TC-P34', 'TD-P35', 'TB-P36', 'TB-P37', 'TC-P38', 'TC-P39', 'TC-P40', 'TC-P41', 'TC-P42', 'TB-P43', 'TD-P44', 'TC-P45', 'TD-P46', 'TD-P47', 'TB-P48', 'TC-P49', 'TD-P50', 'TB-P51', 'TA-P52', 'TA-P53', 'TD-P54', 'TD-P55', 'TB-P56', 'TD-P57', 'TA-P58', 'TC-P59', 'TD-P60', 'TC-P61', 'TC-P62', 'TC-P63', 'TC-P64', 'TD-P65', 'TD-P66', 'TD-P67', 'TD-P68', 'TD-P69', 'TA-P70', 'TB-P71', 'TD-P72', 'TA-P73', 'TD-P74', 'TD-P75', 'TA-P76', 'TC-P77', 'TD-P78', 'TC-P79', 'TC-P80', 'TD-P81', 'TC-P82', 'TD-P83', 'TD-P84', 'TB-P85', 'TC-P86', 'TD-P87', 'TD-P88', 'TD-P89', 'TA-P90', 'TD-P91', 'TD-P92', 'TC-P93', 'TA-P94', 'TC-P95'}, 'RowLabels', table2cell(tasks_name_table(task_significative_change_no_restrictive,:)))%, 'ColumnLabelsRotate', 45)
%print(gcf, 'Task_Analysis_Row_Clustering_2.svg', '-dsvg'); % Save the figure as an SVG file

%% Cluster/Trauma enrichment analysis

clusters = {'cluster_1','cluster_2','cluster_3','cluster_4','cluster_5'};
cluster_trauma_enriched = zeros(4,length(clusters));

for j=1:length(clusters) % 5 = number of clusters
    jth_cluster = eval(clusters{j});
    jth_column_label = jth_cluster.ColumnLabels;
    % Check group A
    count = 0;
    for i=1:length(jth_column_label)
        match = double(find(contains('TA', jth_column_label{i}(1:2))));
        if length(match) > 0
            count =  count + match/length(jth_column_label);
        end
    end
    cluster_trauma_enriched(1,j) = count;
    % Check group B
    count = 0;
    for i=1:length(jth_column_label)
        match = double(find(contains('TB', jth_column_label{i}(1:2))));
        if length(match) > 0
            count =  count + match/length(jth_column_label);
        end    
    end
    cluster_trauma_enriched(2,j) = count;
    % Check group C
    count = 0;
    for i=1:length(jth_column_label)
        match = double(find(contains('TC', jth_column_label{i}(1:2))));
        if length(match) > 0
            count =  count + match/length(jth_column_label);
        end    
    end
    cluster_trauma_enriched(3,j) = count;
    % Check group D
    count = 0;
    for i=1:length(jth_column_label)
        match = double(find(contains('TD', jth_column_label{i}(1:2))));
        if length(match) > 0
            count =  count + match/length(jth_column_label);
        end    
    end
    cluster_trauma_enriched(4,j) = count;
end







cluster_metabogroup_evol = readtable('Task_analysis_Groups_2.xlsx','sheet','Sheet2'); % The ID contains the patient's number and group


% Assuming your table is named 'cluster_metabogroup_evol'
traumaGroups = table2array(cluster_metabogroup_evol(:,2:end));
￼


% Create a figure
figure;

% Plot each trauma group as a separate line
plot(traumaGroups(1,:), 'DisplayName', 'Trauma_group_A'); hold on;
plot(traumaGroups(2,:), 'DisplayName', 'Trauma_group_B'); hold on;
plot(traumaGroups(3,:), 'DisplayName', 'Trauma_group_C'); hold on;
plot(traumaGroups(4,:), 'DisplayName', 'Trauma_group_D'); hold on;

% Add labels and title
xlabel('Cluster');
ylabel('Value');
title('Trauma Groups vs Clusters');

% Add a legend
legend('show');





% Create a figure
figure;

% Plot each trauma group as a separate bar in the bar graph
bar(traumaGroups', 'stacked');

% Add labels and title
xlabel('Cluster');
ylabel('Contribution');
title('Trauma Groups vs Clusters');

% Add a legend
legend(cluster_metabogroup_evol.Properties.RowNames, 'Location', 'northwest');







% Assuming your table is named 'cluster_metabogroup_evol'
traumaGroups = table2array(cluster_metabogroup_evol(:,2:end));

% Create a figure
figure;

% Define a colormap or create a custom color matrix
colormap = [0.9 0 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];

% Plot each trauma group as a separate bar in the bar graph
bar(traumaGroups', 'stacked', 'FaceColor', 'flat');

% Apply the colormap to each bar
for k = 1:size(traumaGroups, 1)
    b(k).CData = colormap(k,:);
end

% Add labels and title
xlabel('Cluster');
ylabel('Value');
title('Trauma Groups vs Clusters');

% Add a legend
legend(cluster_metabogroup_evol.Properties.RowNames, 'Location', 'northwest');



















