%% Load data for plotting MLE results

MLE = load('MLE_result.mat');
sensSoln = load('complex_dusp1_sens.mat').sensSoln;
comp_Model = load('complex_dusp1_model.mat').Model;


%% Making MLE FIM plot with MH points
fims = comp_Model.computeFIM(sensSoln.sens);
[FIM,sFIMcov,fimMetrics] = comp_Model.evaluateExperiment(fims,200);
fimFreePars = FIM(comp_Model.fittingOptions.modelVarsToFit,comp_Model.fittingOptions.modelVarsToFit);
comp_Model.makeMleFimPlot(squeeze(MLE(1,:,:)),fimFreePars,[1,2],0.95)