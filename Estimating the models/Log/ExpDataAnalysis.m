clear 
close all
clc

%% Regression Analysis Script
% In this script there are two significant parts, one part that initialises
% and pre-processes the training data and one part that takes this data and
% feeds it to models for training and testing of machine learning methods.

method =  [true;   % Stepwise Model 1
            false;   %Tree Model 2
            false;   % tree prediction 3
            false;   % PLSR 4
            false;   % GLM 5
            false;   % NLIO 6
            false;   % Optimised Linear Model  7
            false;   % NN  8
            false];    % NARX 9

% saving the models in .mat files        
saveModel = true;

%% Initialize training data
load('datamatrix_log_500')

%pai - annulus pressure, well 1-3
pai_1 = [];  %y1

%pwh - well head pressure, well 1-3
pwh_1 = []; %y4

%wro - wellhead gas production rate, well 1-3
wro_1 = []; %y7

%wrg - wellhead oil production rate, well 1-3
wrg_1 = []; %y10

%prh - riser head pressure
prh = []; %y13

%pm - manifold pressure
pm = []; %y14

%wto - riser head total oil production rate
wto = []; %y15

%wtg - riser head total gas production rate
wtg = [];  %y16

% Erosion rate of change
rocArray = []; %erosion well 1

% Gas Lift Rate
gLift = []; %control input well 1

% Sand Rate
sand = [];

% Time variable 
time = 1:1:501;

% Cumulative Erosion
cumNormErosion = [];

varNames = {'Sand rate'; 'Annulus Pressure'; 'Well head pressure'; 'Well head oil production rate'; 'Well head gas production rate';
    'Riser head pressure';'Manifold pressure';'Riser head total oil production rate';
    'Riser head total gas production rate';'Gas lift rate'};
varNamesShorthand = {'m_s'; 'P_a '; 'P_w_h'; 'WHO'; 'WHG'; 'P_r_h'; 'P_m'; 'RHO'; 'RHG'; 'GLR'; 'ER'};

responseName = {'Erosion rate of change'};

[sandArray,sandArrayNoise,stepSandArray] = sandproductionrate(0.01,501,'log',0.02,50); % The function used for generating logit data
%[sandArray,sandArrayNoise,stepSandArray] = sandproductionrate(0.01,501,'exp',0.005,50); % Function to make exponential data
stepSandArray = modelSandProdRate(stepSandArray,50); %If extrapolated sand rate should be instead of the previous known sand rate

for i = 1:100
    
    % Adding 100 timerseries of yMeas to the end of eachother to give more
    % training data
    
    %yMeas
    pai_1 = [pai_1 Data.yMeas{i,1}(1,:)];
    pwh_1 = [pwh_1 Data.yMeas{i,1}(4,:)];
    wro_1 = [wro_1 Data.yMeas{i,1}(7,:)];
    wrg_1 = [wrg_1 Data.yMeas{i,1}(10,:)];
    prh = [prh Data.yMeas{i,1}(13,:)];
    pm = [pm Data.yMeas{i,1}(14,:)];
    wto = [wto Data.yMeas{i,1}(15,:)];
    wtg = [wtg Data.yMeas{i,1}(16,:)];
    sand = [sand stepSandArray];
    %uArray
    gLift = [gLift Data.uArray{i,1}(1,:)];
        
    %erosionArray - Calculate erosion rate of change
    
    
    erosArray = Data.erosionArray{i,1}(1,:);
    cumNormErosion = [cumNormErosion normalize(erosArray)];
    roc = diff(erosArray)./diff(time);
    roc = [roc, roc(500)];
    rocArray = [rocArray, roc];
    

end

%Data without sand, use if working with constant spr data set
%RegData = [];
%RegData = [transpose(pai_1), transpose(pwh_1), transpose(wro_1), transpose(wrg_1), transpose(prh), transpose(pm), transpose(wto), transpose(wtg), transpose(gLift), transpose(rocArray)];
%NNRegData = [transpose(pai_1), transpose(pwh_1), transpose(wro_1), transpose(wrg_1), transpose(prh), transpose(pm), transpose(wto), transpose(wtg), transpose(gLift)];
%NNRegResponse = [transpose(rocArray)];

%Data with sandproduction rate, use for non-constant spr
RegData = [];
RegData = [transpose(sand), transpose(pai_1), transpose(pwh_1), transpose(wro_1), transpose(wrg_1), transpose(prh), transpose(pm), transpose(wto), transpose(wtg), transpose(gLift), transpose(rocArray)];
NNRegData = [transpose(sand), transpose(pai_1), transpose(pwh_1), transpose(wro_1), transpose(wrg_1), transpose(prh), transpose(pm), transpose(wto), transpose(wtg), transpose(gLift)];
NNRegResponse = [transpose(rocArray)];
NNRegCumResponse = [transpose(cumNormErosion)];
NNRegCumLogarithmicResponse = log(NNRegCumResponse);

regrCenter = mean(NNRegData,1)';
regrScale = std(NNRegData,[],1)';

predCenter = mean(NNRegResponse,1)';
predScale = std(NNRegResponse,[],1)';

NNRegData = normalize(NNRegData);
NNRegResponse = normalize(NNRegResponse);
NormRegData = normalize(RegData);

nameSave = 'normalizationValues';
save(nameSave,'regrCenter','regrScale','predCenter','predScale')

%% Initialize Test Data

%pai - annulus pressure, well 1-3
pai_1_T = [];  %y1

%pwh - well head pressure, well 1-3
pwh_1_T = []; %y4

%wro - wellhead gas production rate, well 1-3
wro_1_T = []; %y7

%wrg - wellhead oil production rate, well 1-3
wrg_1_T = []; %y10

%prh - riser head pressure
prh_T = []; %y13

%pm - manifold pressure
pm_T = []; %y14

%wto - riser head total oil production rate
wto_T = []; %y15

%wtg - riser head total gas production rate
wtg_T = [];  %y16

% Cumulative erosion
cumNormErosionTest = [];

% Erosion rate of change
rocArray_T = []; %erosion well 1
roc_T = [];

% Gas Lift Rate
gLift_T = []; %control input well 1

for i = 101:200
    pai_1_T = [pai_1_T Data.yMeas{i,1}(1,:)];
    pwh_1_T = [pwh_1_T Data.yMeas{i,1}(4,:)];
    wro_1_T = [wro_1_T Data.yMeas{i,1}(7,:)];
    wrg_1_T = [wrg_1_T Data.yMeas{i,1}(10,:)];
    prh_T = [prh_T Data.yMeas{i,1}(13,:)];
    pm_T = [pm_T Data.yMeas{i,1}(14,:)];
    wto_T = [wto_T Data.yMeas{i,1}(15,:)];
    wtg_T = [wtg_T Data.yMeas{i,1}(16,:)];

    gLift_T = [gLift_T Data.uArray{i,1}(1,:)];


 
    erosArray_T = Data.erosionArray{i,1}(1,:);
    cumNormErosionTest = [cumNormErosionTest erosArray_T];
    roc_T = [roc_T gradient(erosArray_T)];

end 
%Use when testing performance of training data with spr
TestData = [transpose(sand), transpose(pai_1_T), transpose(pwh_1_T), transpose(wro_1_T), transpose(wrg_1_T), transpose(prh_T), transpose(pm_T), transpose(wto_T), transpose(wtg_T), transpose(gLift_T)];

%Use when testing performance of training data without spr
%TestData = [transpose(pai_1_T), transpose(pwh_1_T), transpose(wro_1_T), transpose(wrg_1_T), transpose(prh_T), transpose(pm_T), transpose(wto_T), transpose(wtg_T), transpose(gLift_T)];

NormTestResponse = transpose(normalize(roc_T));
NNRegCumResponseTest = [transpose(normalize(cumNormErosionTest))];
NormTestData = normalize(TestData);
NNRegCumLogarithmicResponseTest = log(NNRegCumResponseTest);


%% Training and Testing Stepwise Model
if method(1)
    %Without PCA    
    [StepwiseModel, StepwiseValidationRMSE] = trainStepwiseRegression(NNRegData,NNRegResponse);

    %With PCA
    %[StepwiseModelPCA, StepwisePCAValidationRMSE] = trainStepwiseRegressionPCA(NNRegData,NNRegResponse);

    yhat_stepwise = StepwiseModel.predictFcn(NormTestData);
    %yhat_stepwise_pca = StepwiseModelPCA.predictFcn(NormTestData);

    residuals_stepwise = NormTestResponse - yhat_stepwise;

    mse_stepwise = immse(yhat_stepwise,NormTestResponse)
    %mse_stepwise_pca = immse(yhat_stepwise_pca,NormTestResponse)


    %Plot sample results
    figure(1)

    plot(1:500,yhat_stepwise(1:500),'b')
    hold on
    %plot(1:500,yhat_stepwise_pca(1:500),'b')
    %hold on
    plot(1:500,NormTestResponse(1:500),'r')
    title('Stepwise linear regression')
    legend('Predicted w/ stepwise','True rate')

        %ylim([0,2.2]);
    xlabel('Time [days]');
    ylabel('Erosion Rate (Normalised)');
    
    if saveModel 
        name = 'StepwiseLinear';
        save(name,'StepwiseModel');
    end
end

%% Training and Testing Tree Model
if method(2)
    %Without PCA    
    [OptTreeModel, OptTreeValidationRMSE] = trainOptTreeRegression(NNRegData,NNRegResponse);

    %With PCA
    [OptTreeModelPCA, OptTreePCAValidationRMSE] = trainOptTreeRegressionPCA(NNRegData,NNRegResponse);

    yhat_opttree = OptTreeModel.predictFcn(NormTestData);
    yhat_opttree_pca = OptTreeModelPCA.predictFcn(NormTestData);

    residuals_opttree = NormTestResponse - yhat_opttree;

    mse_opttree = immse(yhat_opttree,NormTestResponse)
    mse_opttree_pca = immse(yhat_opttree_pca,NormTestResponse)

    %Plot sample results
    figure(2)

    plot(1:500,yhat_opttree(1:500),'b')
    hold on
    %plot(1:500,yhat_opttree_pca(1:500),'b')
    %hold on
    plot(1:500,NormTestResponse(1:500),'r')
    title('Optimised Regression Tree')
    legend('Predicted w/ tree','True rate')

        %ylim([0,2.2]);
    xlabel('Time [day]');
    ylabel('Rate of change');
end

% %% Training and Testing SVM Model
% 
% 
% %Without PCA    
% [SvmMedGausModel, SvmMedGausValidationRMSE] = trainSvmMedGausRegression(NNRegData,NNRegResponse);
% 
% %With PCA
% %[SvmMedGausModelPCA, SvmMedGausPCAValidationRMSE] = trainSvmMedGausRegressionPCA(NNRegData,NNRegResponse);
% %%
% %Independent Test
% yhat_svm = optsvm.predictFcn(NormTestData);
% %yhat_svm_pca = SvmMedGausModelPCA.predictFcn(NormTestData);
% 
% residuals_svm = NormTestResponse - yhat_svm;
% 
% %MSE calculation
% mse_svm = immse(yhat_svm,NormTestResponse)
% %mse_svm_pca = immse(yhat_svm_pca,NormTestResponse)


% %% Plot sample results
% figure(3)
%     
% plot(1:500,yhat_svm(1:500),'b')
% hold on
% %plot(1:500,yhat_svm_pca(1:500),'b')
% %hold on
% plot(1:500,NormTestResponse(1:500),'r')
% 
% legend('Predicted w/ Gaussian SVR','True rate')
% title('Medium Gaussian Support Vector Regression')
% %subtitle('with and without PCA')
% 
%     %ylim([0,2.2]);
% xlabel('Time [day]');
% ylabel('Rate of change');

%% Training and Testing Ensemble Model
if method(3)
    %Without PCA
    [EnsembleModel, EnsembleValidationRMSE] = trainEnsembleRegression(NNRegData,NNRegResponse);
    
    %With PCA
    [EnsembleModelPCA, EnsemblePCAValidationRMSE] = trainEnsembleRegressionPCA(NNRegData,NNRegResponse);
    
    %Independent Test
    yhat_ensemble = EnsembleModel.predictFcn(NormTestData);
    yhat_ensemble_pca = EnsembleModelPCA.predictFcn(NormTestData);
    
    residuals_ensemble = NormTestResponse - yhat_ensemble;
    
    %MSE calculation
    mse_ensemble = immse(yhat_ensemble,NormTestResponse)
    mse_ensemble_pca = immse(yhat_ensemble_pca,NormTestResponse)
    
    %Plot sample results
    yhat_ensemble_low = yhat_ensemble - 0.102*1.65;
    yhat_ensemble_high = yhat_ensemble + 0.102*1.65;
    
    figure(4)
    plot(1:500,yhat_ensemble(1:500),'b')
    hold on
    %plot(1:500,yhat_ensemble_pca(1:500),'b')
    %hold on
    plot(1:500,NormTestResponse(1:500),'r')
    %hold on
    %plot(1:500,yhat_ensemble_low(1:500),'g',1:500,yhat_ensemble_high(1:500),'g')
    title('Optimised ensemble using modelled sand production rate')
    legend('Bagged ensemble prediction','True rate')
    
    %ylim([0,2.2]);
    xlabel('Time [day]');
    ylabel('Rate of change');
    
    %Plotting residuals and normdist
    figure(12)
    yyaxis left
    title('Residuals of Ensemble Predictions and Normal PDF')
    %subtitle('SD=0.102, mean = 0.0069')
    histogram(residuals_ensemble);
    hold on
    yyaxis right
    x_values = -0.5:0.01:0.5;
    y = pdf(ensembleDist,x_values);
    plot(x_values,y,'r')
end

%% PLSR 
if method(4)
    X = NormRegData(:,1:10);
    y = NormRegData(:,11);   

    [XL,yl,XS,YS,beta,PCTVAR,mse,stats] = plsregress(X,y,5);

    %plot(1:9,cumsum(100*PCTVAR(2,:)),'-bo');
    %xlabel('Number of PLS components');
    %ylabel('Percent Variance Explained in y');

    yfit = [ones(size(X,1),1) X]*beta;
    residuals_plsr = y - yfit;
    %stem(residuals)

    Xtestplot = NormTestData(1:500,1:10);
    yfitTestPlot = [ones(size(Xtestplot,1),1) Xtestplot]*beta;

    mse_plsr = immse(yfit,NormTestResponse)

    figure(5)

    plot(1:500,yfitTestPlot(1:500),'b')
    hold on
    plot(1:500,NormTestResponse(1:500),'r')
    title('Partial Least Squares')
    legend('PLSR','True rate')

        %ylim([0,2.2]);
    xlabel('Time [day]');
    ylabel('Rate of change');


    figure(6)
    plot(1:5 ,cumsum(100*PCTVAR(2,:)),'-bo');
    xlabel('Number of PLS components');
    ylabel('Percent Variance Explained in Y');


    figure(7)
    plot(1:10,stats.W(:,1:2),'-');
    xlabel('Variable');
    ylabel('PLS Weight');
    legend({'1st Component' '2nd Component' '3rd Component'},'location','NW');
end
%% GLM Modelling
if method(5)
    [b,dev,stats] = glmfit(NNRegData,NNRegResponse,'normal');
    yhat = glmval(b,NormTestData,'identity');

    mse_glm = immse(yhat,NormTestResponse)


    figure(8)
    plot(1:500,NormTestResponse(1:500),'r',1:500,yhat(1:500),'b')
    title('GLM predictions and True values')
    legend({'Measured','GLM Prediction'})

end
%% Neural Network Modelling (Non-linear Input/Output)
if method(6)
    NLIO_Logdata;


    x1 = transpose(NormTestData(4:50100,:));
    xi1 = transpose(NormTestData(1:3,:));


    [y1,xf1] = NLIOnet(x1,xi1)

    mse_nlio = immse(transpose(y1),NormTestResponse(4:50100))

    figure(9)
    plot(4:500,NormTestResponse(4:500),'r',4:500,y1(4:500),'b')
    title('Non-linear I/O Neural Network')
    legend({'Measured','Neural Network Prediction'})
end

%% Optimised Linear Model
if method(7)
    hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus');
    [Mdl,FitInfo,HyperparameterOptimizationResults] = fitrlinear(NNRegData,NNRegResponse,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',hyperopts)

    [yfit_regularised] = predict(Mdl,NormTestData);

    mse_ridge = immse(yfit_regularised,NormTestResponse)

    figure(10)
    plot(1:500,NormTestResponse(1:500),'r',1:500,yfit_regularised(1:500),'b')
    title('Optimised Regularised Linear Regression')
    %subtitle('Least squares, expected-improvement, ridge regularisation')
    legend({'Measured','Regularised Linear Regression'})
end

%% Neural Network Fitting
if method(8)
    NNfit_nonTS;

    yfit_nnfit = net(transpose(NormTestData));
    mse_nnfit = immse(transpose(yfit_nnfit),NormTestResponse)

    figure(10)
    plot(1:500,NormTestResponse(1:500),'r',1:500,yfit_nnfit(1:500),'b')
    title('Artificial Neural Network Regression')
    xlabel('Time [days]')
    ylabel('Erosion Rate (norm.)')
    legend({'Measured','ANN regression fit'})
    
   %if saveModel 
        %name = 'ANN';
        %For NN use script: NNFit.m
    %end
    
end

%% Experimenting with cumulative erosion modelling with NARX
if method(9)
    NARX_LogExpDataCumulative;
    
    X = tonndata(NormTestData(1:10000,:),false,false);
    T = tonndata(NNRegCumResponseTest(1:10000),false,false);
    
    [x,xi,ai,t] = preparets(net,X,{},T);
    
    numTimesteps = size(x,2);
    knownOutputTimesteps = 1:(numTimesteps-250);
    predictOutputTimesteps = (numTimesteps-249):numTimesteps;
    X1 = X(:,knownOutputTimesteps);
    T1 = T(:,knownOutputTimesteps);
    [x1,xio,aio] = preparets(net,X1,{},T1);
    [y1,xfo,afo] = net(x1,xio,aio);
    
    x2 = X(1,predictOutputTimesteps);
    [netc,xic,aic] = closeloop(net,xfo,afo);
    [y2,xfc,afc] = netc(x2,xic,aic);
    
    y1 = cell2mat(y1);
    y2 = cell2mat(y2);
    
    mse_narx = immse(transpose(y2),NNRegCumResponseTest(9751:10000))
    
    figure(6)
    plot(9556:9750,transpose(y1(9550:9744)),'-o')
    hold on
    plot(9751:10000,transpose(y2),'-o')
    hold on
    plot(9550:10000,NNRegCumResponseTest(9550:10000),'r')
    title('NARX Neural Network Prediction and true values')
    %subtitle('Cumulative error at endpoint 0.0392')
    legend({'Open Loop Prediction','Closed Loop Prediction','True Values'})
end

%% Function Definitions of model training functions

% stepwise functions w and w/o pca
function [trainedModel, validationRMSE] = trainStepwiseRegression(trainingData, responseData)


inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

concatenatedPredictorsAndResponse = predictors;
concatenatedPredictorsAndResponse.NNRegResponse = response;
linearModel = stepwiselm(...
    concatenatedPredictorsAndResponse, ...
    'linear', ...
    'Upper', 'interactions', ...
    'NSteps', 1000, ...
    'Verbose', 0);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
linearModelPredictFcn = @(x) predict(linearModel, x);
trainedModel.predictFcn = @(x) linearModelPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedModel.LinearModel = linearModel;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 10 columns because this model was trained using 10 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
KFolds = 5;
cvp = cvpartition(size(response, 1), 'KFold', KFolds);
% Initialize the predictions to the proper sizes
validationPredictions = response;
for fold = 1:KFolds
    trainingPredictors = predictors(cvp.training(fold), :);
    trainingResponse = response(cvp.training(fold), :);
    foldIsCategoricalPredictor = isCategoricalPredictor;
    
    % Train a regression model
    % This code specifies all the model options and trains the model.
    concatenatedPredictorsAndResponse = trainingPredictors;
    concatenatedPredictorsAndResponse.NNRegResponse = trainingResponse;
    linearModel = stepwiselm(...
        concatenatedPredictorsAndResponse, ...
        'linear', ...
        'Upper', 'interactions', ...
        'NSteps', 1000, ...
        'Verbose', 0);
    
    % Create the result struct with predict function
    linearModelPredictFcn = @(x) predict(linearModel, x);
    validationPredictFcn = @(x) linearModelPredictFcn(x);
    
    % Add additional fields to the result struct
    
    % Compute validation predictions
    validationPredictors = predictors(cvp.test(fold), :);
    foldPredictions = validationPredictFcn(validationPredictors);
    
    % Store predictions in the original order
    validationPredictions(cvp.test(fold), :) = foldPredictions;
end

% Compute validation RMSE
isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));

end

function [trainedModel, validationRMSE] = trainStepwiseRegressionPCA(trainingData, responseData)

inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
isCategoricalPredictorBeforePCA = isCategoricalPredictor;
numericPredictors = predictors(:, ~isCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% 'inf' values have to be treated as missing data for PCA.
numericPredictors(isinf(numericPredictors)) = NaN;
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = [array2table(pcaScores(:,1:numComponentsToKeep)), predictors(:, isCategoricalPredictor)];
isCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(isCategoricalPredictor))];

% Train a regression model
% This code specifies all the model options and trains the model.
concatenatedPredictorsAndResponse = predictors;
concatenatedPredictorsAndResponse.NNRegResponse = response;
linearModel = stepwiselm(...
    concatenatedPredictorsAndResponse, ...
    'linear', ...
    'Upper', 'interactions', ...
    'NSteps', 1000, ...
    'Verbose', 0);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
linearModelPredictFcn = @(x) predict(linearModel, x);
trainedModel.predictFcn = @(x) linearModelPredictFcn(pcaTransformationFcn(predictorExtractionFcn(x)));

% Add additional fields to the result struct
trainedModel.PCACenters = pcaCenters;
trainedModel.PCACoefficients = pcaCoefficients;
trainedModel.LinearModel = linearModel;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 10 columns because this model was trained using 10 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
KFolds = 5;
cvp = cvpartition(size(response, 1), 'KFold', KFolds);
% Initialize the predictions to the proper sizes
validationPredictions = response;
for fold = 1:KFolds
    trainingPredictors = predictors(cvp.training(fold), :);
    trainingResponse = response(cvp.training(fold), :);
    foldIsCategoricalPredictor = isCategoricalPredictor;
    
    % Apply a PCA to the predictor matrix.
    % Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
    isCategoricalPredictorBeforePCA = foldIsCategoricalPredictor;
    numericPredictors = trainingPredictors(:, ~foldIsCategoricalPredictor);
    numericPredictors = table2array(varfun(@double, numericPredictors));
    % 'inf' values have to be treated as missing data for PCA.
    numericPredictors(isinf(numericPredictors)) = NaN;
    [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
        numericPredictors);
    % Keep enough components to explain the desired amount of variance.
    explainedVarianceToKeepAsFraction = 95/100;
    numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
    pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
    trainingPredictors = [array2table(pcaScores(:,1:numComponentsToKeep)), trainingPredictors(:, foldIsCategoricalPredictor)];
    foldIsCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(foldIsCategoricalPredictor))];
    
    % Train a regression model
    % This code specifies all the model options and trains the model.
    concatenatedPredictorsAndResponse = trainingPredictors;
    concatenatedPredictorsAndResponse.NNRegResponse = trainingResponse;
    linearModel = stepwiselm(...
        concatenatedPredictorsAndResponse, ...
        'linear', ...
        'Upper', 'interactions', ...
        'NSteps', 1000, ...
        'Verbose', 0);
    
    % Create the result struct with predict function
    pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
    linearModelPredictFcn = @(x) predict(linearModel, x);
    validationPredictFcn = @(x) linearModelPredictFcn(pcaTransformationFcn(x));
    
    % Add additional fields to the result struct
    
    % Compute validation predictions
    validationPredictors = predictors(cvp.test(fold), :);
    foldPredictions = validationPredictFcn(validationPredictors);
    
    % Store predictions in the original order
    validationPredictions(cvp.test(fold), :) = foldPredictions;
end

% Compute validation RMSE
isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));
end

% optimised tree functions w and w/o pca
function [trainedModel, validationRMSE] = trainOptTreeRegression(trainingData, responseData)

inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Train a regression model
% This code specifies all the model options and trains the model.
regressionTree = fitrtree(...
    predictors, ...
    response, ...
    'MinLeafSize', 68, ...
    'Surrogate', 'off');

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
treePredictFcn = @(x) predict(regressionTree, x);
trainedModel.predictFcn = @(x) treePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedModel.RegressionTree = regressionTree;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 10 columns because this model was trained using 10 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedModel.RegressionTree, 'KFold', 5);

% Compute validation predictions
validationPredictions = kfoldPredict(partitionedModel);

% Compute validation RMSE
validationRMSE = sqrt(kfoldLoss(partitionedModel, 'LossFun', 'mse'));
end

function [trainedModel, validationRMSE] = trainOptTreeRegressionPCA(trainingData, responseData)

inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
isCategoricalPredictorBeforePCA = isCategoricalPredictor;
numericPredictors = predictors(:, ~isCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% 'inf' values have to be treated as missing data for PCA.
numericPredictors(isinf(numericPredictors)) = NaN;
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = [array2table(pcaScores(:,1:numComponentsToKeep)), predictors(:, isCategoricalPredictor)];
isCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(isCategoricalPredictor))];

% Train a regression model
% This code specifies all the model options and trains the model.
regressionTree = fitrtree(...
    predictors, ...
    response, ...
    'MinLeafSize', 28, ...
    'Surrogate', 'off');

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
treePredictFcn = @(x) predict(regressionTree, x);
trainedModel.predictFcn = @(x) treePredictFcn(pcaTransformationFcn(predictorExtractionFcn(x)));

% Add additional fields to the result struct
trainedModel.PCACenters = pcaCenters;
trainedModel.PCACoefficients = pcaCoefficients;
trainedModel.RegressionTree = regressionTree;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 10 columns because this model was trained using 10 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
KFolds = 5;
cvp = cvpartition(size(response, 1), 'KFold', KFolds);
% Initialize the predictions to the proper sizes
validationPredictions = response;
for fold = 1:KFolds
    trainingPredictors = predictors(cvp.training(fold), :);
    trainingResponse = response(cvp.training(fold), :);
    foldIsCategoricalPredictor = isCategoricalPredictor;
    
    % Apply a PCA to the predictor matrix.
    % Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
    isCategoricalPredictorBeforePCA = foldIsCategoricalPredictor;
    numericPredictors = trainingPredictors(:, ~foldIsCategoricalPredictor);
    numericPredictors = table2array(varfun(@double, numericPredictors));
    % 'inf' values have to be treated as missing data for PCA.
    numericPredictors(isinf(numericPredictors)) = NaN;
    [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
        numericPredictors);
    % Keep enough components to explain the desired amount of variance.
    explainedVarianceToKeepAsFraction = 95/100;
    numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
    pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
    trainingPredictors = [array2table(pcaScores(:,1:numComponentsToKeep)), trainingPredictors(:, foldIsCategoricalPredictor)];
    foldIsCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(foldIsCategoricalPredictor))];
    
    % Train a regression model
    % This code specifies all the model options and trains the model.
    regressionTree = fitrtree(...
        trainingPredictors, ...
        trainingResponse, ...
        'MinLeafSize', 28, ...
        'Surrogate', 'off');
    
    % Create the result struct with predict function
    pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
    treePredictFcn = @(x) predict(regressionTree, x);
    validationPredictFcn = @(x) treePredictFcn(pcaTransformationFcn(x));
    
    % Add additional fields to the result struct
    
    % Compute validation predictions
    validationPredictors = predictors(cvp.test(fold), :);
    foldPredictions = validationPredictFcn(validationPredictors);
    
    % Store predictions in the original order
    validationPredictions(cvp.test(fold), :) = foldPredictions;
end

% Compute validation RMSE
isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));
end

% medium gaussian svm w and w/o pca
function [trainedModel, validationRMSE] = trainSvmMedGausRegression(trainingData, responseData)

inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Train a regression model
% This code specifies all the model options and trains the model.
responseScale = iqr(response);
if ~isfinite(responseScale) || responseScale == 0.0
    responseScale = 1.0;
end
boxConstraint = responseScale/1.349;
epsilon = responseScale/13.49;
regressionSVM = fitrsvm(...
    predictors, ...
    response, ...
    'KernelFunction', 'gaussian', ...
    'PolynomialOrder', [], ...
    'KernelScale', 3.2, ...
    'BoxConstraint', boxConstraint, ...
    'Epsilon', epsilon, ...
    'Standardize', true);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
svmPredictFcn = @(x) predict(regressionSVM, x);
trainedModel.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedModel.RegressionSVM = regressionSVM;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 10 columns because this model was trained using 10 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
KFolds = 5;
cvp = cvpartition(size(response, 1), 'KFold', KFolds);
% Initialize the predictions to the proper sizes
validationPredictions = response;
for fold = 1:KFolds
    trainingPredictors = predictors(cvp.training(fold), :);
    trainingResponse = response(cvp.training(fold), :);
    foldIsCategoricalPredictor = isCategoricalPredictor;
    
    % Train a regression model
    % This code specifies all the model options and trains the model.
    responseScale = iqr(trainingResponse);
    if ~isfinite(responseScale) || responseScale == 0.0
        responseScale = 1.0;
    end
    boxConstraint = responseScale/1.349;
    epsilon = responseScale/13.49;
    regressionSVM = fitrsvm(...
        trainingPredictors, ...
        trainingResponse, ...
        'KernelFunction', 'gaussian', ...
        'PolynomialOrder', [], ...
        'KernelScale', 3.2, ...
        'BoxConstraint', boxConstraint, ...
        'Epsilon', epsilon, ...
        'Standardize', true);
    
    % Create the result struct with predict function
    svmPredictFcn = @(x) predict(regressionSVM, x);
    validationPredictFcn = @(x) svmPredictFcn(x);
    
    % Add additional fields to the result struct
    
    % Compute validation predictions
    validationPredictors = predictors(cvp.test(fold), :);
    foldPredictions = validationPredictFcn(validationPredictors);
    
    % Store predictions in the original order
    validationPredictions(cvp.test(fold), :) = foldPredictions;
end

% Compute validation RMSE
isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));
end

function [trainedModel, validationRMSE] = trainSvmMedGausRegressionPCA(trainingData, responseData)

inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
isCategoricalPredictorBeforePCA = isCategoricalPredictor;
numericPredictors = predictors(:, ~isCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% 'inf' values have to be treated as missing data for PCA.
numericPredictors(isinf(numericPredictors)) = NaN;
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = [array2table(pcaScores(:,1:numComponentsToKeep)), predictors(:, isCategoricalPredictor)];
isCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(isCategoricalPredictor))];

% Train a regression model
% This code specifies all the model options and trains the model.
responseScale = iqr(response);
if ~isfinite(responseScale) || responseScale == 0.0
    responseScale = 1.0;
end
boxConstraint = responseScale/1.349;
epsilon = responseScale/13.49;
regressionSVM = fitrsvm(...
    predictors, ...
    response, ...
    'KernelFunction', 'gaussian', ...
    'PolynomialOrder', [], ...
    'KernelScale', 3.2, ...
    'BoxConstraint', boxConstraint, ...
    'Epsilon', epsilon, ...
    'Standardize', true);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
svmPredictFcn = @(x) predict(regressionSVM, x);
trainedModel.predictFcn = @(x) svmPredictFcn(pcaTransformationFcn(predictorExtractionFcn(x)));

% Add additional fields to the result struct
trainedModel.PCACenters = pcaCenters;
trainedModel.PCACoefficients = pcaCoefficients;
trainedModel.RegressionSVM = regressionSVM;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 10 columns because this model was trained using 10 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
KFolds = 5;
cvp = cvpartition(size(response, 1), 'KFold', KFolds);
% Initialize the predictions to the proper sizes
validationPredictions = response;
for fold = 1:KFolds
    trainingPredictors = predictors(cvp.training(fold), :);
    trainingResponse = response(cvp.training(fold), :);
    foldIsCategoricalPredictor = isCategoricalPredictor;
    
    % Apply a PCA to the predictor matrix.
    % Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
    isCategoricalPredictorBeforePCA = foldIsCategoricalPredictor;
    numericPredictors = trainingPredictors(:, ~foldIsCategoricalPredictor);
    numericPredictors = table2array(varfun(@double, numericPredictors));
    % 'inf' values have to be treated as missing data for PCA.
    numericPredictors(isinf(numericPredictors)) = NaN;
    [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
        numericPredictors);
    % Keep enough components to explain the desired amount of variance.
    explainedVarianceToKeepAsFraction = 95/100;
    numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
    pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
    trainingPredictors = [array2table(pcaScores(:,1:numComponentsToKeep)), trainingPredictors(:, foldIsCategoricalPredictor)];
    foldIsCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(foldIsCategoricalPredictor))];
    
    % Train a regression model
    % This code specifies all the model options and trains the model.
    responseScale = iqr(trainingResponse);
    if ~isfinite(responseScale) || responseScale == 0.0
        responseScale = 1.0;
    end
    boxConstraint = responseScale/1.349;
    epsilon = responseScale/13.49;
    regressionSVM = fitrsvm(...
        trainingPredictors, ...
        trainingResponse, ...
        'KernelFunction', 'gaussian', ...
        'PolynomialOrder', [], ...
        'KernelScale', 3.2, ...
        'BoxConstraint', boxConstraint, ...
        'Epsilon', epsilon, ...
        'Standardize', true);
    
    % Create the result struct with predict function
    pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
    svmPredictFcn = @(x) predict(regressionSVM, x);
    validationPredictFcn = @(x) svmPredictFcn(pcaTransformationFcn(x));
    
    % Add additional fields to the result struct
    
    % Compute validation predictions
    validationPredictors = predictors(cvp.test(fold), :);
    foldPredictions = validationPredictFcn(validationPredictors);
    
    % Store predictions in the original order
    validationPredictions(cvp.test(fold), :) = foldPredictions;
end

% Compute validation RMSE
isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));
end

% optimised ensemble functions w and w/o pca
function [trainedModel, validationRMSE] = trainEnsembleRegression(trainingData, responseData)

% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Train a regression model
% This code specifies all the model options and trains the model.
template = templateTree(...
    'MinLeafSize', 4, ...
    'NumVariablesToSample', 8);
regressionEnsemble = fitrensemble(...
    predictors, ...
    response, ...
    'Method', 'Bag', ...
    'NumLearningCycles', 408, ...
    'Learners', template);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
ensemblePredictFcn = @(x) predict(regressionEnsemble, x);
trainedModel.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedModel.RegressionEnsemble = regressionEnsemble;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 10 columns because this model was trained using 10 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedModel.RegressionEnsemble, 'KFold', 5);

% Compute validation predictions
validationPredictions = kfoldPredict(partitionedModel);

% Compute validation RMSE
validationRMSE = sqrt(kfoldLoss(partitionedModel, 'LossFun', 'mse'));
end

function [trainedModel, validationRMSE] = trainEnsembleRegressionPCA(trainingData, responseData)

inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Apply a PCA to the predictor matrix.
% Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
isCategoricalPredictorBeforePCA = isCategoricalPredictor;
numericPredictors = predictors(:, ~isCategoricalPredictor);
numericPredictors = table2array(varfun(@double, numericPredictors));
% 'inf' values have to be treated as missing data for PCA.
numericPredictors(isinf(numericPredictors)) = NaN;
[pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
    numericPredictors);
% Keep enough components to explain the desired amount of variance.
explainedVarianceToKeepAsFraction = 95/100;
numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
predictors = [array2table(pcaScores(:,1:numComponentsToKeep)), predictors(:, isCategoricalPredictor)];
isCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(isCategoricalPredictor))];

% Train a regression model
% This code specifies all the model options and trains the model.
template = templateTree(...
    'MinLeafSize', 772, ...
    'NumVariablesToSample', 6);
regressionEnsemble = fitrensemble(...
    predictors, ...
    response, ...
    'Method', 'LSBoost', ...
    'NumLearningCycles', 456, ...
    'Learners', template, ...
    'LearnRate', 0.9564018983943808);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
ensemblePredictFcn = @(x) predict(regressionEnsemble, x);
trainedModel.predictFcn = @(x) ensemblePredictFcn(pcaTransformationFcn(predictorExtractionFcn(x)));

% Add additional fields to the result struct
trainedModel.PCACenters = pcaCenters;
trainedModel.PCACoefficients = pcaCoefficients;
trainedModel.RegressionEnsemble = regressionEnsemble;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 10 columns because this model was trained using 10 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9', 'column_10'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false];

% Perform cross-validation
KFolds = 5;
cvp = cvpartition(size(response, 1), 'KFold', KFolds);
% Initialize the predictions to the proper sizes
validationPredictions = response;
for fold = 1:KFolds
    trainingPredictors = predictors(cvp.training(fold), :);
    trainingResponse = response(cvp.training(fold), :);
    foldIsCategoricalPredictor = isCategoricalPredictor;
    
    % Apply a PCA to the predictor matrix.
    % Run PCA on numeric predictors only. Categorical predictors are passed through PCA untouched.
    isCategoricalPredictorBeforePCA = foldIsCategoricalPredictor;
    numericPredictors = trainingPredictors(:, ~foldIsCategoricalPredictor);
    numericPredictors = table2array(varfun(@double, numericPredictors));
    % 'inf' values have to be treated as missing data for PCA.
    numericPredictors(isinf(numericPredictors)) = NaN;
    [pcaCoefficients, pcaScores, ~, ~, explained, pcaCenters] = pca(...
        numericPredictors);
    % Keep enough components to explain the desired amount of variance.
    explainedVarianceToKeepAsFraction = 95/100;
    numComponentsToKeep = find(cumsum(explained)/sum(explained) >= explainedVarianceToKeepAsFraction, 1);
    pcaCoefficients = pcaCoefficients(:,1:numComponentsToKeep);
    trainingPredictors = [array2table(pcaScores(:,1:numComponentsToKeep)), trainingPredictors(:, foldIsCategoricalPredictor)];
    foldIsCategoricalPredictor = [false(1,numComponentsToKeep), true(1,sum(foldIsCategoricalPredictor))];
    
    % Train a regression model
    % This code specifies all the model options and trains the model.
    template = templateTree(...
        'MinLeafSize', 772, ...
        'NumVariablesToSample', 6);
    regressionEnsemble = fitrensemble(...
        trainingPredictors, ...
        trainingResponse, ...
        'Method', 'LSBoost', ...
        'NumLearningCycles', 456, ...
        'Learners', template, ...
        'LearnRate', 0.9564018983943808);
    
    % Create the result struct with predict function
    pcaTransformationFcn = @(x) [ array2table((table2array(varfun(@double, x(:, ~isCategoricalPredictorBeforePCA))) - pcaCenters) * pcaCoefficients), x(:,isCategoricalPredictorBeforePCA) ];
    ensemblePredictFcn = @(x) predict(regressionEnsemble, x);
    validationPredictFcn = @(x) ensemblePredictFcn(pcaTransformationFcn(x));
    
    % Add additional fields to the result struct
    
    % Compute validation predictions
    validationPredictors = predictors(cvp.test(fold), :);
    foldPredictions = validationPredictFcn(validationPredictors);
    
    % Store predictions in the original order
    validationPredictions(cvp.test(fold), :) = foldPredictions;
end

% Compute validation RMSE
isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));
end