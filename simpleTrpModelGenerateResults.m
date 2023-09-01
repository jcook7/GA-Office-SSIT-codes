%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
  clear all
  clc
%% Define SSIT Model
  Model = SSIT;
  Model.species = {'x1';'x2'}; % x1:gene state, x2:mRNA
  Model.initialCondition = [0;0];
  Model.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
  Model.stoichiometry = [1,-1,0,0;0,0,1,-1];
  Model.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))'};
  Model.parameters = ({'koff',0.14;'kon',0.14;'kr',25;'gr',0.01;...
                       'a1',0.4;'r1',0.04;'r2',0.1});
  Model.fspOptions.initApproxSS = true;
  %% Tryptolide Experiment
ModelTrypt = Model;
ModelTrypt.propensityFunctions(3) = {'kr*x1*Itrypt'};
ModelTrypt.inputExpressions(2,:) = {'Itrypt','(t<tpt)'};
tpt_array = 20:20:180;
ModelTrypt.tSpan = tpt_array;
ModelTrypt.sensOptions.solutionMethod = 'finiteDifference';

ModelTrypt.solutionScheme = 'FSP';
ModelTrypt.fspOptions.fspTol = 1e-6;
ModelTrypt.parameters(8,:) = {'tpt',180};
ModelTrypt.fspOptions.bounds=[];
[fspSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;
ModelTrypt.fspOptions.bounds

%%  Setting data
ModelTrypt = Model;
% load dataset
ModelTrypt = ModelTrypt.loadData('../ExampleData/DUSP1_Dex_100nM_Rep1_Rep2.csv',{'x3','RNA_nuc'});
ModelTrypt.fittingOptions.modelVarsToFit = 1:8;

ModelTrypt.propensityFunctions(3) = {'kr*x1*Itrypt'};
ModelTrypt.inputExpressions(2,:) = {'Itrypt','(t<tpt)'};
tpt_array = 20:20:180;
ModelTrypt.sensOptions.solutionMethod = 'finiteDifference';


%% Running FSP fits
for i=1:5
    ModelTrypt.solutionScheme = 'FSP';
    ModelTrypt.fspOptions.fspTol = 1e-6;
    ModelTrypt.parameters(8,:) = {'tpt',180};
    ModelTrypt.fspOptions.bounds=[];
    [fspSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;
    ModelTrypt.fspOptions.bounds

% Load and Fit smFISH Data
      ModelTrypt.fspOptions.fspTol = inf;
      ModelTrypt.fittingOptions.modelVarsToFit = 1:8;
      fitOptions = optimset('Display','iter','MaxIter',500);
      ModelTrypt.parameters(1:8,2) = num2cell(ModelTrypt.maximizeLikelihood([],fitOptions));
      ModelTrypt.makeFitPlot;
end


%% Metropolis Hastings to Quantify Parameter Uncertainty (TRP)
  ModelTrypt.fittingOptions.modelVarsToFit = 1:8;
  for i =1:5
      MHOptions = struct('numberOfSamples',1000,'burnin',0,'thin',1,...
          'useFIMforMetHast',true,'suppressFSPExpansion',true);
      [bestParsFound,~,mhResults] = ModelTrypt.maximizeLikelihood([ModelTrypt.parameters{ModelTrypt.fittingOptions.modelVarsToFit,2}]',...
          MHOptions,'MetropolisHastings');
      ModelTrypt.parameters(ModelTrypt.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
      ModelTrypt.plotMHResults(mhResults);
  end
%% Save results
[sensSoln,ModelTrypt.fspOptions.bounds] = ModelTrypt.solve;