%% Full Model (EGRNT)
tmp = load('complex_dusp1_model');
EGRNT = tmp.Model;

EGRNT.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
EGRNT.fspOptions.fspTol = 1e-8;  % Set FSP error tolerance.
[EGRNTsoln,EGRNT.fspOptions.bounds] = EGRNT.solve;  % Solve the FSP analysis
[EGRNTsoln,EGRNT.fspOptions.bounds] = EGRNT.solve;  % Solve the FSP analysis

% EGRNT.ssaOptions.nSimsPerExpt = 100;
% EGRNT.ssaOptions.Nexp = 200; 
% EGRNT.sampleDataFromFSP(EGRNTsoln,'full_dusp1_model_testC.csv'); 

EGRNT.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
%EGRNT = EGRNT.loadData('full_dusp1_model_testC.csv',{'x3','exp1_s3'});
EGRNT = EGRNT.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x3','RNA_nuc'});
EGRNT.initialTime = 0;
EGRNT.fittingOptions.timesToFit = ones(1,length(EGRNT.tSpan),'logical');
EGRNT.makeFitPlot

tic
for i=1:5
[EGRNTsoln,EGRNT.fspOptions.bounds] = EGRNT.solve;  % Solve the FSP analysis
end
complexComputeTime = toc;
complexComputeTime/5

%% Reduced model (SGRS)
clc
SGRS = SSIT;
SGRS.species = {'x1';'x2'};  % GRnuc, geneOn, dusp1
SGRS.initialCondition = [0;0];
SGRS.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
SGRS.inputExpressions = {'IGR','kcn0/knc+(t>=0)*kcn1/(r1-knc)*(exp(-knc*t)-exp(-r1*t))'};
SGRS.stoichiometry = [1,-1,0,0;0,0,1,-1];
SGRS.parameters = EGRNT.parameters;
SGRS.tSpan = EGRNT.tSpan;
SGRS.fspOptions.initApproxSS = true;

tmp2 = load('simple_dusp1_model');
SGRS = tmp2.simple_Model;
SGRS.parameters=EGRNT.parameters;

SGRS.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
SGRS.fspOptions.fspTol = 1e-8;  % Set FSP error tolerance.
[SGRSsoln,SGRS.fspOptions.bounds] = SGRS.solve;  % Solve the FSP analysis
[SGRSsoln,SGRS.fspOptions.bounds] = SGRS.solve;  % Solve the FSP analysis
SGRS.computeLikelihood()

tic
for i=1:5
[SGRSsoln,SGRS.fspOptions.bounds] = SGRS.solve;  % Solve the FSP analysis
end
simpleComputeTime = toc;
simpleComputeTime/5

%% Plot Comparison of Full and Reduce Model
SGRS.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
%SGRS = SGRS.loadData('full_dusp1_model_testC.csv',{'x2','exp1_s3'});
SGRS = SGRS.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x3','RNA_nuc'});

SGRS.initialTime = 0;
SGRS.fittingOptions.timesToFit = ones(1,length(SGRS.tSpan),'logical');
SGRS.makeFitPlot

% to get the likelihood for the loaded data:
SGRS.computeLikelihood

%%
clc, clear
simple_Model = load('simple_dusp1_model.mat').simple_Model;
cellCounts = simple_Model.dataSet.nCells;
simple_sensSoln = load('simple_dusp1_sens.mat').sensSoln;
comp_Model = load('complex_dusp1_model.mat').Model;
comp_sensSoln = load('complex_dusp1_sens.mat').sensSoln;
comp_fimResults = comp_Model.computeFIM(comp_sensSoln.sens);
%%

simplefimResults = simple_Model.computeFIM(simple_sensSoln.sens);
simple_Model.pdoOptions.unobservedSpecies = {'x1'};
[simpleFIM,sFIMcov1,fimMetrics1] = simple_Model.evaluateExperiment(simplefimResults,cellCounts);
iFIMs=inv(simpleFIM);


comp_sensSoln = load('complex_dusp1_sens.mat').sensSoln;
comp_Model.pdoOptions.unobservedSpecies = {'x1','x2'};
[compFIM,sFIMcov2,fimMetrics2] = comp_Model.evaluateExperiment(comp_fimResults,cellCounts);
iFIMc =inv(compFIM);


varSimple = var(diag(iFIMs))
varComplex = var(diag(iFIMc))
