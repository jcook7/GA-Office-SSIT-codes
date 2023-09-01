%% Practice Challenge Script
% care most about the gene in the on state
% 1.) It is optimal to identify what species are in the diagram and from there make a map of all possible reactions we care about in our system
% 2.) Then make propensities
% 3.) Finally do stiocheometry
% kon is correlated to the gene being in the on state
% koff is correlated to the gene being in the off state
% kr is associated with the rate of transcription of gene 1 (species 1) for Nascent mRNA (species 2)
% gene 1 (species 1) is in the on-state which allows Nascent mRNA (species 2) to be transcribed from it (mRNA being produced at a rate (kr))
% 'kr = k0 + (k1 / (alpha + (D^eta)))' where:
%       k0 is the unbound state
%       k1 is the bound state
%       alpha is the initiation rate
%       D is the burst duration
%       eta is likely the burst frequency
% gamma is the degredation rate of Nascent mRNA (species 2)
Example = SSIT();    % Create SSIT instance using pre-selected model
Example.species = {'x1';'x2'};   % Set species names. [Exercise Note: x1 == gene 1 (species 1), x2 == species 2]
Example.parameters = {'kon',5;'koff',10;'k0',1;'k1',1;'alpha',1;'D',1;'eta',1;'gamma',1}; % Set parameter names and values
Example.propensityFunctions = {'kon*(2-x1)';'koff*(x1)';...
    'x1*(k0 + (k1 / (alpha + (D^eta))))';'gamma*x2'}; % Set propensity functions
% In this case species 'x1' is the number of alleles in the 'on' state. There is a maximum of two alleles
% so the propensity function of the gene state switching to an active allele is "kon*(2-x1)".
% When the system has two alleles in the 'on' state, the propensity function of gaining another 'on' state
% is zero due to a maximum of two possible alleles being in the on state.
% What reaction can occur at an infinitesimally small time span? With a lot of mRNA present, it is likely that a degradation reaction will occur.
% Therefore, gamma*(x2 being species indicating the amount of mRNA present) is the degredation propensity function.
Example.stoichiometry = [1,-1,0, 0;
                         0, 0,1,-1]; % Columns are reactions and number of rows are the numbers of associated species
                                     % Row 1 is species 1 rxns
                                     % Row 2 is species 2 rxns
                                     % The first rxn for gene 1 (species 1) will have to be associated with the gene turning on (indicated by the positive 1 in the first column of the first row)
                                     % The second rxn for gene 1 (species 1) will have to be associated with the gene turning off (indicated by the negative 1 in the second column of the first row)
                                     % gene 1 (species 1) is in the on-state which allows Nascent mRNA (species 2) to be transcribed from it (mRNA being produced at a rate (kr)) (indicated by the positive 1 in the third column of the second row)
                                     % gamma is associated with the degredation rate of Nascent mRNA (species 2) (indicated by the negative 1 in the fourth column of the second row)
Example.initialCondition = [0;0]; % Set initial condition
% Example.customConstraintFuns = {'(x1-3).^2.*(x2-3).^2'}; % Set FSP constraint. We don't know anything about is
Example.tSpan = [1000:50:2000];        % Set times at which to compute distributions (look at time span in CSV file to determine correct timespan)
%% Loading Data
Example = Example.loadData('../ExampleData/drug3_TitrationB.csv',{'x2','num_rna_nuc'});
% Because we care about fitting our nascent RNA (which is species x2 in this case) and the associated Nascent mRNA in the nucleus in our csv file is called
% 'num_rna_nuc', so we need to input the species and its associated column
% You can check this by putting 'Example' in the command window. This helps
% you navigate your model.
% Example.dataSet will show you how the data is loaded into the model
% If anything goes wrong you may not have any cells loaded in or the link
% species and data may not make any sense (like a 0x0 [rows by columns])
% other things to look at are the times. Something went wrong if the times
% don't match up to the csv files.
%% Solve using the FSP approach
Example.solutionScheme = 'FSP';    % Set solutions scheme to FSP. [Tells us we want to explicitly solve for the FPS in this case.]
Example.fspOptions.fspTol = 1e-4;  % Set FSP error tolerance. [Associated with the truncation. Error tolerance basically states we know what the expected error of using the actual chemical master equation vs. using the FSP is going to be, so we can create a more accurate model]
[FSPsoln,Example.fspOptions.bounds] = Example.solve;  % Solve the FSP analysis
Example.makePlot(FSPsoln,'marginals',[2:5],false,[1,2])    % Plot marginal distributions
Example.makePlot(FSPsoln,'joints',[2:5],false,[5])         % Plot joint distributions

%% FSP fitting
for i=1:10
    Example.fittingOptions.modelVarsToFit = 1:8
    Example.solutionScheme = 'FSP';
    Example.fspOptions.fspTol = 1e-4;
    Example.fspOptions.verbose = 0;
    Example.fspOptions.bounds=[];
    [fspSoln,Example.fspOptions.bounds] = Example.solve;
    Example.fspOptions.bounds
    
    % Load and Fit smFISH Data
    Example.fspOptions.fspTol = inf;
    fitOptions = optimset('Display','iter','MaxIter',500);
    Example.parameters(Example.fittingOptions.modelVarsToFit,2) = num2cell(Example.maximizeLikelihood([Example.parameters{Example.fittingOptions.modelVarsToFit,2}],fitOptions));
end
Example.makeFitPlot;

%% Metropolis Hastings to Quantify Parameter Uncertainty
Example.fittingOptions.modelVarsToFit = 1:9;
MHOptions = struct('numberOfSamples',1500,'burnin',100,'thin',1,...
  'useFIMforMetHast',true,'suppressFSPExpansion',true);
allFitOptions.CovFIMscale = 0.1;
[bestParsFound,~,mhResults] = Example.maximizeLikelihood([Example.parameters{Example.fittingOptions.modelVarsToFit,2}]',...
  MHOptions,'MetropolisHastings');
Example.parameters(Example.fittingOptions.modelVarsToFit,2) = num2cell(bestParsFound);
Example.plotMHResults(mhResults);

%% Calculate Sensitivity Mtarix
Example.sensOptions.solutionMethod = 'finiteDifference';
Example.solutionScheme = 'fspSens';
Example.fspOptions.fspTol = 1e-6;
[sensSoln,Example.fspOptions.bounds] = Example.solve;

%% Calculate FIM
fims = Example.computeFIM(sensSoln.sens);
numberCellsPerTimePoint = Example.dataSet.nCells;
FIM = Example.evaluateExperiment(fims,numberCellsPerTimePoint);
