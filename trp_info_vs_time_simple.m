%% Using the SSIT to fit and Design Single-Cell Experiments 
% In this script, we show how the SSIT can be used to identify a
% time-inhomogeneous model for the activation of Dusp1 mRNA expression
% under Dexamethasome stimulation of Glucocorticoid Receptors.
  clear all
  clc
%% Simple SSIT Model
  Model = SSIT;
  Model.species = {'x1';'x2'}; % x1:gene state, x2:mRNA
  Model.initialCondition = [0;0];
  Model.propensityFunctions = {'kon*IGR*(2-x1)';'koff*x1';'kr*x1';'gr*x2'};
  Model.stoichiometry = [1,-1,0,0;0,0,1,-1];
  Model.inputExpressions = {'IGR','1+a1*exp(-r1*t)*(1-exp(-r2*t))'};
  Model.parameters = ({'koff',0.14;'kon',0.14;'kr',25;'gr',0.01;...
                       'a1',0.4;'r1',0.04;'r2',0.1});
  Model.fspOptions.initApproxSS = true;

  %% Load saved mHast and Sensitivity

  mhResults = load('simple_dusp1_Trp_mhast.mat').mhResults;
  sensSoln = load('simple_dusp1_Trp_sens.mat').sensSoln;
  ModelTrypt = load('simple_dusp1_Trp_model.mat').ModelTrypt;

%% |COV| plotted vs time
  %trpt_times_vect = linspace(0,180,10);
  trpt_times_vect = [0:10:50,55:5:180];
  %trpt_times_vect = [0:10:40,42:2:180];
  for i = length(trpt_times_vect):-1:1
      i
      tic
      ModelTrypt_t{i} = ModelTrypt;
      ModelTrypt_t{i}.parameters{9,2} = trpt_times_vect(i);
      ModelTrypt_t{i}.sensOptions.solutionMethod = 'finiteDifference';
      ModelTrypt_t{i}.solutionScheme = 'fspSens';
      ModelTrypt_t{i}.fspOptions.fspTol = 1e-6;
      [sensSoln_t{i},ModelTrypt_t{i}.fspOptions.bounds] = ModelTrypt_t{i}.solve;
      fimResults{i} = ModelTrypt_t{i}.computeFIM(sensSoln_t{i}.sens);
      FIM_all{i} = ModelTrypt_t{i}.evaluateExperiment(fimResults{i},ModelTrypt_t{i}.dataSet.nCells);
      FIM_det(i) = det(FIM_all{i}(1:7,1:7));
      
      toc
  end 
  
  %% Plot
  figure(1)
  plot(trpt_times_vect,FIM_det.^(-1))
  title('|COV| vs time')
  ylabel('|COV|')
  xlabel('measurement timepoints [min]')


  %% Optimized fims
  for i = length(fimResults):-1:1
      nTotal = sum(ModelTrypt_t{i}.dataSet.nCells);
      nCellsOpt = ModelTrypt_t{i}.optimizeCellCounts(fimResults{i},nTotal,'TR[1:4]');
      fimOpt{i} = ModelTrypt_t{i}.evaluateExperiment(fimResults{i},nCellsOpt);
      FIM_det_opt(i) = det(FIM_all{i}(1:8,1:8));
  end
    
  figure(2)
  plot(trpt_times_vect,FIM_det_opt.^(-1))
  title('|COV| vs time [Optimized]')
  ylabel('|COV|')
  xlabel('measurement timepoints [min]')
