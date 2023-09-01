%% example_BPPB -- example script to show SSIT for the BPPB Seminar
clear all
close all
clc
%% Complex Dusp1 model
Model = SSIT;
Model.species = {'x1';'x2';'x3'};  % GRnuc, geneOn, dusp1
Model.initialCondition = [0;0;0];
Model.propensityFunctions = {'(kcn0+kcn1*IDex)';'knc*x1';...
  'kon*x1*(2-x2)';'koff*x2';'kr*x2';'gr*x3'};
Model.inputExpressions = {'IDex','(t>0)*exp(-r1*t)'};
Model.parameters = ({'koff',0.14;'kon',0.01;'kr',1;'gr',0.08;...
  'kcn0',0.01;'kcn1',0.1;'knc',1;'r1',0.01});
Model.stoichiometry = [ 1,-1, 0, 0, 0, 0;...
                      0, 0, 1,-1, 0, 0;...
                      0, 0, 0, 0, 1,-1];
Model.fspOptions.initApproxSS = true;


A = load('complex_dusp1_model.mat').Model;
%A = load('simple_dusp1_model.mat').simple_Model;
%A.fittingOptions.logPrior=[];

%% simple model
% Define SSIT Model
  Model = SSIT;
  Model.species = {'x1';'x2'};
  Model.initialCondition = [0;0];
  Model.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
  Model.stoichiometry = [1,-1,0,0;0,0,1,-1];
  Model.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))'};
  Model.parameters = ({'koff',0.14;'kon',0.14;'kr',5;'gr',0.5;...
                       'a1',0.4;'r1',0.04;'r2',0.1});
  Model.fspOptions.initApproxSS = true;

  A = load('simple_dusp1_model.mat').simple_Model;

%% Create SSIT Model
% A = SSIT();    % Create SSIT instance using pre-selected model
% A.species = {'x1';'x2'};   % Set species names.
% A.parameters = {'kb',10;'ka',80;'M',20;'g',1}; % Set parameter names and values
% A.stoichiometry = [1,-1,0, 0; 0, 0,1,-1]; % Set Stoichiometry matrix
% A.propensityFunctions = {'kb+ka*M^3/(M^3+x2^3)';'g*x1';...
%     'kb+ka*M^3/(M^3+x1^3)';'g*x2'}; % Set propensity functions
% A.initialCondition = [4;0]; % Set initial condition
% A.customConstraintFuns = {'(x1-3).^2.*(x2-3).^2'}; % Set FSP constraint.
% A.tSpan = [0:3:12];          % Set times at which to compute distributions

%% Solve using the FSP approach
A.solutionScheme = 'FSP';    % Set solutions scheme to FSP.
A.fspOptions.fspTol = 1e-4;  % Set FSP error tolerance.
[FSPsoln,A.fspOptions.bounds] = A.solve;  % Solve the FSP analysis
%A.makePlot(FSPsoln,'marginals',[2:5],false,[1,2])    % Plot marginal distributions
%A.makePlot(FSPsoln,'joints',[2:5],false,[5])         % Plot joint distributions

%% Solve Sensitivity using FSP
A.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln] = A.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem
%A.makePlot(sensSoln,'marginals',[],false,[3,4]) % Plot marginal sensitivities

%% Define a Binomial PDO
%A.pdoOptions.type = 'Binomial';
%A.pdoOptions.props.CaptureProbabilityS1 = 0.7;  % Distortion for S1
%A.pdoOptions.props.CaptureProbabilityS2 = 0.7; % Distortion for S2
%A.pdoOptions.PDO = A.generatePDO(A.pdoOptions,[0.7,0.7],FSPsoln.fsp,true);
%figure(20); contourf(log10(A.pdoOptions.PDO.conditionalPmfs{1})); colorbar
%xlabel('"true" number of mRNA'); ylabel('observed number of mRNA'); set(gca,'fontsize',15);

%% Apply PDO to FSP and Sensitivity Calculations
A.solutionScheme = 'FSP'; % Set solution scheme to FSP.
%A.makePlot(FSPsoln,'marginals',[2:5],true,[1,2])  % Plot Distorted Marginals
A.solutionScheme = 'fspSens'; % Set solution scheme to Sensitivity
%A.makePlot(sensSoln,'marginals',[],true,[3,4])    % Plot Distorted Sensitivities

%% Simulate Data for Subsequent Fitting
comp_Modelssa=A;
comp_Modelssa.solutionScheme = 'SSA';  % Set solution scheme to SSA.
comp_Modelssa.ssaOptions.Nexp = 10; comp_Modelssa.ssaOptions.nSimsPerExpt = 200;
%A.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
A.solve([],'ToggleSSAData50Expts.csv'); 
%% 
starttime=-15000
comp_Modelssa=A;
comp_Modelssa.initialTime = starttime;
comp_Modelssa.tSpan = [starttime,A.tSpan];
comp_Modelssa.ssaOptions.useTimeVar=true;
comp_Modelssa.ssaOptions.signalUpdateRate=1;
comp_Modelssa.solutionScheme = 'SSA';  % Set solution scheme to SSA.
comp_Modelssa.ssaOptions.Nexp = 10; 
comp_Modelssa.ssaOptions.nSimsPerExpt = 200
A.ssaOptions.applyPDO = true; % Include the distortion in the SSA data.
comp_Modelssa.solve([],'complex_model_test.csv') 
comp_Modelssa.solutionScheme = 'FSP'  % Set solution scheme back to FSP.
%% Associate Datasets with FSP Models.
A.solutionScheme = 'FSP';  % Set solution scheme back to FSP.
%B = A.loadData('simple_model_test.csv',{'x1','exp1_s1';'x2','exp1_s2'});
B = A.loadData('simple_model_test.csv',{'x1','exp1_s1';'x2','exp1_s2';'x3','exp1_s3'});

%B{1}.pdoOptions.PDO = [];  % Do not use PDO.

%% Plot Results
%for m=1:3; [fit_error,fitSolutuions] = B{m}.computeLikelihood; B{m}.makeFitPlot(fitSolutuions); end

%% Find MLE for each simulated data set.
B.fittingOptions.modelVarsToFit = [1:4];
nfreepar = length(B.fittingOptions.modelVarsToFit);
MLE = zeros(1,nfreepar,comp_Modelssa.ssaOptions.Nexp); 
fMLE = inf(1,nfreepar,comp_Modelssa.ssaOptions.Nexp);
for iExp = 1:comp_Modelssa.ssaOptions.Nexp
    iExp
    %B = B.loadData('simple_model_test.csv',{'x1',['exp',num2str(iExp),'_s1'];'x2',['exp',num2str(iExp),'_s2']}); % Link non-distorted data.
    B = B.loadData('simple_model_test.csv',{'x1',['exp',num2str(iExp),'_s1'];'x2',['exp',num2str(iExp),'_s2'];'x3',['exp',num2str(iExp),'_s3']}); % Link non-distorted data.
    %B = B.loadData('ToggleSSAData50Expts.csv',{'x1',['exp',num2str(iExp),'_s1']}); % Link non-distorted data.000
    B.fittingOptions.timesToFit = [false,ones(1,length(A.tSpan),'logical')];
    B.tSpan = B.tSpan(2:end);
    % Link non-distorted data.
        
        if iExp==1; x0 = [B.parameters{B.fittingOptions.modelVarsToFit,2}]';
        else; x0 = squeeze(MLE(1,:,iExp-1)); end
        fitOptions = optimset('Display','iter','MaxIter',10);
        [MLE(1,:,iExp),fMLE(1,:,iExp)] = B.maximizeLikelihood(x0,fitOptions);
    
end

%% Compute FIM
cellCounts = B.ssaOptions.nSimsPerExpt*ones(size(B.tSpan));  % Number of cells in each experiment.
B.pdoOptions.unobservedSpecies = {};
fimResults = B.computeFIM(sensSoln.sens); 
% Compute the FIM for full observations and no distortion.
[FIM,sFIMcov,fimMetrics] = B.evaluateExperiment(fimResults,cellCounts);
fimFreePars = FIM(B.fittingOptions.modelVarsToFit,B.fittingOptions.modelVarsToFit);
B.makeMleFimPlot(squeeze(MLE(1,:,1:10)),fimFreePars,[2,4],0.95)


%% Optimize Experiment
A.tSpan = [0:1:99]; % Set times at which to compute distributions
A.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
[sensSoln,bounds] = A.solve;  % Solve the sensitivity problem
cellCounts = 100*ones(size(A.tSpan));
nCellsTotal = sum(cellCounts);
% Compute FIM with all parameters free to vary
A.fittingOptions.pdoVarsToFit = 'all';
fimResults = A.computeFIM(sensSoln.sens); 
FIMintuit = A.evaluateExperiment(fimResults,cellCounts);
% Optimize cell counts for different objectives
nCellsOptDet = A.optimizeCellCounts(fimResults,nCellsTotal,'Determinant');
FIMoptDet= A.evaluateExperiment(fimResults,nCellsOptDet);
nCellsOpt1234 = A.optimizeCellCounts(fimResults,nCellsTotal,'[1:4]');
FIMoptDet1234= A.evaluateExperiment(fimResults,nCellsOpt1234);

% Plot results
pars = [[A.parameters{:,2}],[0.7,0.7]]';
for i=1:5; for j = i+1:6
        subplot(5,5,(i-1)*5+j-1)
        A.makeMleFimPlot([],FIMintuit,[j,i],0.95,1,pars([j,i])); hold on
        A.makeMleFimPlot([],FIMoptDet,[j,i],0.95,1,pars([j,i]))
        A.makeMleFimPlot([],FIMoptDet1234,[j,i],0.95,1,pars([j,i]))
        legend("off")
end;end
legend('Intuitive','D-Optimality','$\min \det(\Sigma_{\rm model})$','interpreter','latex')
