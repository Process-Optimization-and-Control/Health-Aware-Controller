%% Regression Analysis Script
% In this script there are two significant parts, one part that initialises
% and pre-processes the training data and one part that takes this data and
% feeds it to models for training and testing of machine learning methods.

%% Initialize training data

clc 

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

% Cumulative erosion

cumNormErosion = [];

varNames = {'Sand rate'; 'Annulus Pressure'; 'Well head pressure'; 'Well head oil production rate'; 'Well head gas production rate';
    'Riser head pressure';'Manifold pressure';'Riser head total oil production rate';
    'Riser head total gas production rate';'Gas lift rate'};
varNamesShorthand = {'P_a '; 'P_w_h'; 'WHO'; 'WHG'; 'P_r_h'; 'P_m'; 'RHO'; 'RHG'; 'GLR'; 'ER'};

responseName = {'Erosion rate of change'};

%[sandArray,sandArrayNoise,stepSandArray] = sandproductionrate(0.01,501,'exp',0.005);

%% 
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
    %sand = [sand stepSandArray];
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
RegData = [];
RegData = [transpose(pai_1), transpose(pwh_1), transpose(wro_1), transpose(wrg_1), transpose(prh), transpose(pm), transpose(wto), transpose(wtg), transpose(gLift), transpose(rocArray)];
NNRegData = [transpose(pai_1), transpose(pwh_1), transpose(wro_1), transpose(wrg_1), transpose(prh), transpose(pm), transpose(wto), transpose(wtg), transpose(gLift)];
NNRegResponse = [transpose(rocArray)];
NNRegCumResponse = [transpose(cumNormErosion)];


%Data with sandproduction rate, use for non-constant spr
%RegData = [];
%RegData = [transpose(sand), transpose(pai_1), transpose(pwh_1), transpose(wro_1), transpose(wrg_1), transpose(prh), transpose(pm), transpose(wto), transpose(wtg), transpose(gLift), transpose(rocArray)];
%NNRegData = [transpose(sand), transpose(pai_1), transpose(pwh_1), transpose(wro_1), transpose(wrg_1), transpose(prh), transpose(pm), transpose(wto), transpose(wtg), transpose(gLift)];
%NNRegResponse = [transpose(rocArray)];

NNRegData = normalize(NNRegData);
NNRegResponse = normalize(NNRegResponse);
NormRegData = normalize(RegData);

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

% Erosion rate of change

rocArray_T = []; %erosion well 1
roc_T = [];

% Gas Lift Rate

gLift_T = []; %control input well 1

% Cumulative Erosion

cumNormErosionTest = [];

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
    roc_T = [roc_T gradient(erosArray_T)];
    cumNormErosionTest = [cumNormErosionTest normalize(erosArray_T)];

end 
%Use when testing performance of training data with spr
TestData = [transpose(sand), transpose(pai_1_T), transpose(pwh_1_T), transpose(wro_1_T), transpose(wrg_1_T), transpose(prh_T), transpose(pm_T), transpose(wto_T), transpose(wtg_T), transpose(gLift_T)];

%Use when testing performance of training data without spr
%TestData = [transpose(pai_1_T), transpose(pwh_1_T), transpose(wro_1_T), transpose(wrg_1_T), transpose(prh_T), transpose(pm_T), transpose(wto_T), transpose(wtg_T), transpose(gLift_T)];

NormTestResponse = transpose(normalize(roc_T));
    
NormTestData = normalize(TestData);
NNRegCumResponseTest = [transpose(cumNormErosionTest)];


%% Training and Testing Stepwise Model

%Without PCA    
[StepwiseModel, StepwiseValidationRMSE] = trainStepwiseRegression(NNRegData,NNRegResponse);

%With PCA
[StepwiseModelPCA, StepwisePCAValidationRMSE] = trainStepwiseRegressionPCA(NNRegData,NNRegResponse);

yhat_stepwise = StepwiseModel.predictFcn(NormTestData);
yhat_stepwise_pca = StepwiseModelPCA.predictFcn(NormTestData);


mse_stepwise = immse(yhat_stepwise,NormTestResponse)
mse_stepwise_pca = immse(yhat_stepwise_pca,NormTestResponse)


%% Plot sample results
figure(1)
    
plot(1:500,yhat_stepwise(1:500),'g')
hold on
plot(1:500,yhat_stepwise_pca(1:500),'b')
hold on
plot(1:500,NormTestResponse(1:500),'r')

legend('Predicted w/ stepwise','Predicted w/ stepwise+pca','True rate')
    
    %ylim([0,2.2]);
xlabel('Time [day]');
ylabel('Rate of change');


% %% Training and Testing Tree Model
% 
% %Without PCA    
% [OptTreeModel, OptTreeValidationRMSE] = trainOptTreeRegression(NNRegData,NNRegResponse);
% 
% %With PCA
% [OptTreeModelPCA, OptTreePCAValidationRMSE] = trainOptTreeRegressionPCA(NNRegData,NNRegResponse);
% 
% yhat_opttree = OptTreeModel.predictFcn(NormTestData);
% yhat_opttree_pca = OptTreeModelPCA.predictFcn(NormTestData);
% 
% 
% mse_opttree = immse(yhat_opttree,NormTestResponse);
% mse_opttree_pca = immse(yhat_opttree_pca,NormTestResponse);
% 
% % Plot sample results
% figure(2)
%     
% plot(1:500,yhat_opttree(1:500),'g')
% hold on
% plot(1:500,yhat_opttree_pca(1:500),'b')
% hold on
% plot(1:500,NormTestResponse(1:500),'r')
% 
% legend('Predicted w/ tree','Predicted w/ tree+pca','True rate')
%     
%     %ylim([0,2.2]);
% xlabel('Time [day]');
% ylabel('Rate of change');

% %% Training and Testing SVM Model
% 
% 
% %Without PCA    
% [SvmMedGausModel, SvmMedGausValidationRMSE] = trainSvmMedGausRegression(NNRegData,NNRegResponse);
% 
% %With PCA
% [SvmMedGausModelPCA, SvmMedGausPCAValidationRMSE] = trainSvmMedGausRegressionPCA(NNRegData,NNRegResponse);
% 
% %Independent Test
% yhat_svm = SvmMedGausModel.predictFcn(NormTestData);
% yhat_svm_pca = SvmMedGausModelPCA.predictFcn(NormTestData);
% 
% %MSE calculation
% mse_svm = immse(yhat_svm,NormTestResponse)
% mse_svm_pca = immse(yhat_svm_pca,NormTestResponse)
% 
% 
% % Plot sample results
% figure(3)
% plot(1:500,yhat_svm(1:500),'g')
% hold on
% plot(1:500,yhat_svm_pca(1:500),'b')
% hold on
% plot(1:500,NormTestResponse(1:500),'r')
% 
% legend('Predicted w/ Gaussian SVR','Predicted w/ Gaussian SVR+PCA','True rate')
%     
%     %ylim([0,2.2]);
% xlabel('Time [day]');
% ylabel('Rate of change');

% %% Training and Testing Ensemble Model
% 
% 
% %Without PCA    
% [EnsembleModel, EnsembleValidationRMSE] = trainEnsembleRegression(NNRegData,NNRegResponse);
% 
% %With PCA
% [EnsembleModelPCA, EnsemblePCAValidationRMSE] = trainEnsembleRegressionPCA(NNRegData,NNRegResponse);
% 
% %Independent Test
% yhat_ensemble = EnsembleModel.predictFcn(NormTestData);
% yhat_ensemble_pca = EnsembleModelPCA.predictFcn(NormTestData);
% 
% %MSE calculation
% mse_ensemble = immse(yhat_ensemble,NormTestResponse)
% mse_ensemble_pca = immse(yhat_ensemble_pca,NormTestResponse)
% 
% % Plot sample results
% figure(4)
%     
% plot(1:500,yhat_ensemble(1:500),'g')
% hold on
% plot(1:500,yhat_ensemble_pca(1:500),'b')
% hold on
% plot(1:500,NormTestResponse(1:500),'r')
% 
% legend('Predicted w/ ensemble','Predicted w/ ensemble+pca','True rate')
%     
%     %ylim([0,2.2]);
% xlabel('Time [day]');
% ylabel('Rate of change');
% 
% %% PLSR 
% 
% X = NormRegData(:,1:9);
% y = NormRegData(:,10);   
% 
% [XL,yl,XS,YS,beta,PCTVAR] = plsregress(X,y,5);
% 
% %plot(1:9,cumsum(100*PCTVAR(2,:)),'-bo');
% %xlabel('Number of PLS components');
% %ylabel('Percent Variance Explained in y');
% 
% yfit = [ones(size(X,1),1) X]*beta;
% residuals = y - yfit;
% %stem(residuals)
% 
% 
% Xtest = NormTestData(:,1:9);
% yfitTest = [ones(size(Xtest,1),1) Xtest]*beta;
% 
% mse_plsr = immse(yfitTest,NormTestResponse)
% 
% figure(5)
%     
% plot(1:500,yfitTest(1:500),'b')
% hold on
% plot(1:500,NormTestResponse(1:500),'r')
% title('Partial Least Squares')
% legend('PLSR','True rate')
%     
%     %ylim([0,2.2]);
% xlabel('Time [day]');
% ylabel('Rate of change');

% %% GLM Modelling
% 
% [b,dev,stats] = glmfit(NNRegData,NNRegResponse,'normal');
% yhat = glmval(b,NormTestData,'identity');
% 
% mse_glm = immse(yhat,NormTestResponse)
% 
% 
% figure(6)
% plot(1:500,NormTestResponse(1:500),'r',1:500,yhat(1:500),'b')
% title('GLM predictions and True values')
% legend({'Measured','GLM Prediction'})

% 
% %% Neural Network Modelling (Non-linear Input/Output)
% 
% NLIO_Logdata;
% 
% x1 = transpose(NormTestData(4:50100,:));
% xi1 = transpose(NormTestData(1:3,:));
% 
% 
% [y1,xf1] = NLIOnetConst(x1,xi1);
% 
% mse_nlio = immse(transpose(y1),NormTestResponse(4:50100))
% 
% figure(7)
% plot(504:1000,NormTestResponse(504:1000),'r',504:1000,transpose(y1(504:1000)),'b')
% title('Non-linear I/O Neural Network')
% legend({'Measured','Neural Network Prediction'})

% %% ANN Regression
% 
% NNFit_constData;
% 
% yfit_nnfit = net(transpose(NormTestData));
% mse_nnfit = immse(transpose(yfit_nnfit),NormTestResponse)
% 
% figure(10)
% plot(1:500,NormTestResponse(1:500),'r',1:500,yfit_nnfit(1:500),'b')
% title('Artificial Neural Network Regression')
% legend({'Measured','ANN regression fit'})

% %% Optimised Linear Model
% hyperopts = struct('AcquisitionFunctionName','expected-improvement-plus');
% [Mdl,FitInfo,HyperparameterOptimizationResults] = fitrlinear(NNRegData,NNRegResponse,'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',hyperopts)
% 
% 
% yfit_regularised = predict(Mdl,NormTestData);
% 
% mse_ridge = immse(yfit_regularised,NormTestResponse)
% 
% figure(8)
% plot(1:500,NormTestResponse(1:500),'r',1:500,yfit_regularised(1:500),'b')
% title('Optimised Regularised Linear Regression')
% subtitle('Least squares, expected-improvement, ridge regularisation')
% legend({'Measured','Regularised Linear Regression'})



% %% Experimenting with cumulative erosion modelling with NARX
% 
% NARX_LogExpDataCumulative; 
% 
% X = tonndata(NormTestData(1:10000,:),false,false);
% T = tonndata(NNRegCumResponseTest(1:10000),false,false);
% 
% [x,xi,ai,t] = preparets(net,X,{},T);
% 
% numTimesteps = size(x,2);
% knownOutputTimesteps = 1:(numTimesteps-250);
% predictOutputTimesteps = (numTimesteps-249):numTimesteps;
% X1 = X(:,knownOutputTimesteps);
% T1 = T(:,knownOutputTimesteps);
% [x1,xio,aio] = preparets(net,X1,{},T1);
% [y1,xfo,afo] = net(x1,xio,aio);
% 
% x2 = X(1,predictOutputTimesteps);
% [netc,xic,aic] = closeloop(net,xfo,afo);
% [y2,xfc,afc] = netc(x2,xic,aic);
% 
% y1 = cell2mat(y1);
% y2 = cell2mat(y2);
% 
% mse_narx = immse(transpose(y2),NNRegCumResponseTest(9751:10000))
% 
% figure(6)
% plot(9556:9750,transpose(y1(9550:9744)),'-o')
% hold on
% plot(9751:10000,transpose(y2),'-o')
% hold on
% plot(9550:10000,NNRegCumResponseTest(9550:10000),'r')
% title('NARX Neural Network Prediction and true values')
% subtitle('Mean Square Error of closed loop prediction: 0.0138')
% legend({'Open Loop Prediction','Closed Loop Prediction','True Values'})


%% Function Definitions of model training functions

% stepwise functions w and w/o pca
function [trainedModel, validationRMSE] = trainStepwiseRegression(trainingData, responseData)
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% Returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: A matrix with the same number of columns and data type
%       as the matrix imported into the app.
%
%      responseData: A vector with the same data type as the vector
%       imported into the app. The length of responseData and the number of
%       rows of trainingData must be equal.
%
%  Output:
%      trainedModel: A struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: A function to make predictions on new data.
%
%      validationRMSE: A double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input arguments trainingData and responseData.
%
% For example, to retrain a regression model trained with the original data
% set T and response Y, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

% Auto-generated by MATLAB on 10-Nov-2020 16:35:58


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

% Train a regression model
% This code specifies all the model options and trains the model.
concatenatedPredictorsAndResponse = predictors;
concatenatedPredictorsAndResponse.ConstNNRegResponse = response;
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
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 9 columns because this model was trained using 9 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
    concatenatedPredictorsAndResponse.ConstNNRegResponse = trainingResponse;
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
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% Returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: A matrix with the same number of columns and data type
%       as the matrix imported into the app.
%
%      responseData: A vector with the same data type as the vector
%       imported into the app. The length of responseData and the number of
%       rows of trainingData must be equal.
%
%  Output:
%      trainedModel: A struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: A function to make predictions on new data.
%
%      validationRMSE: A double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input arguments trainingData and responseData.
%
% For example, to retrain a regression model trained with the original data
% set T and response Y, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

% Auto-generated by MATLAB on 10-Nov-2020 16:41:48


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
concatenatedPredictorsAndResponse.ConstNNRegResponse = response;
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
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 9 columns because this model was trained using 9 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
    concatenatedPredictorsAndResponse.ConstNNRegResponse = trainingResponse;
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
    
    % Compute validation RMSE
isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));
end
end

% optimised tree functions w and w/o pca
function [trainedModel, validationRMSE] = trainOptTreeRegression(trainingData, responseData)
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% Returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: A matrix with the same number of columns and data type
%       as the matrix imported into the app.
%
%      responseData: A vector with the same data type as the vector
%       imported into the app. The length of responseData and the number of
%       rows of trainingData must be equal.
%
%  Output:
%      trainedModel: A struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: A function to make predictions on new data.
%
%      validationRMSE: A double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input arguments trainingData and responseData.
%
% For example, to retrain a regression model trained with the original data
% set T and response Y, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

% Auto-generated by MATLAB on 10-Nov-2020 16:55:42


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

% Train a regression model
% This code specifies all the model options and trains the model.
regressionTree = fitrtree(...
    predictors, ...
    response, ...
    'MinLeafSize', 33, ...
    'Surrogate', 'off');

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
treePredictFcn = @(x) predict(regressionTree, x);
trainedModel.predictFcn = @(x) treePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedModel.RegressionTree = regressionTree;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 9 columns because this model was trained using 9 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedModel.RegressionTree, 'KFold', 5);

% Compute validation predictions
validationPredictions = kfoldPredict(partitionedModel);

% Compute validation RMSE
validationRMSE = sqrt(kfoldLoss(partitionedModel, 'LossFun', 'mse'));
end

function [trainedModel, validationRMSE] = trainOptTreeRegressionPCA(trainingData, responseData)
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% Returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: A matrix with the same number of columns and data type
%       as the matrix imported into the app.
%
%      responseData: A vector with the same data type as the vector
%       imported into the app. The length of responseData and the number of
%       rows of trainingData must be equal.
%
%  Output:
%      trainedModel: A struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: A function to make predictions on new data.
%
%      validationRMSE: A double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input arguments trainingData and responseData.
%
% For example, to retrain a regression model trained with the original data
% set T and response Y, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

% Auto-generated by MATLAB on 10-Nov-2020 16:58:40


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
    'MinLeafSize', 17, ...
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
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 9 columns because this model was trained using 9 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
        'MinLeafSize', 17, ...
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
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% Returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: A matrix with the same number of columns and data type
%       as the matrix imported into the app.
%
%      responseData: A vector with the same data type as the vector
%       imported into the app. The length of responseData and the number of
%       rows of trainingData must be equal.
%
%  Output:
%      trainedModel: A struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: A function to make predictions on new data.
%
%      validationRMSE: A double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input arguments trainingData and responseData.
%
% For example, to retrain a regression model trained with the original data
% set T and response Y, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

% Auto-generated by MATLAB on 10-Nov-2020 17:12:24


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
    'KernelScale', 3, ...
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
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 9 columns because this model was trained using 9 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
        'KernelScale', 3, ...
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
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% Returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: A matrix with the same number of columns and data type
%       as the matrix imported into the app.
%
%      responseData: A vector with the same data type as the vector
%       imported into the app. The length of responseData and the number of
%       rows of trainingData must be equal.
%
%  Output:
%      trainedModel: A struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: A function to make predictions on new data.
%
%      validationRMSE: A double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input arguments trainingData and responseData.
%
% For example, to retrain a regression model trained with the original data
% set T and response Y, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

% Auto-generated by MATLAB on 10-Nov-2020 17:13:43


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
    'KernelScale', 3, ...
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
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 9 columns because this model was trained using 9 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
        'KernelScale', 3, ...
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
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% Returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: A matrix with the same number of columns and data type
%       as the matrix imported into the app.
%
%      responseData: A vector with the same data type as the vector
%       imported into the app. The length of responseData and the number of
%       rows of trainingData must be equal.
%
%  Output:
%      trainedModel: A struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: A function to make predictions on new data.
%
%      validationRMSE: A double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input arguments trainingData and responseData.
%
% For example, to retrain a regression model trained with the original data
% set T and response Y, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

% Auto-generated by MATLAB on 10-Nov-2020 17:17:50


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

% Train a regression model
% This code specifies all the model options and trains the model.
template = templateTree(...
    'MinLeafSize', 17, ...
    'NumVariablesToSample', 9);
regressionEnsemble = fitrensemble(...
    predictors, ...
    response, ...
    'Method', 'LSBoost', ...
    'NumLearningCycles', 29, ...
    'Learners', template, ...
    'LearnRate', 0.3449936310235667);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
ensemblePredictFcn = @(x) predict(regressionEnsemble, x);
trainedModel.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedModel.RegressionEnsemble = regressionEnsemble;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2020b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 9 columns because this model was trained using 9 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

% Perform cross-validation
partitionedModel = crossval(trainedModel.RegressionEnsemble, 'KFold', 5);

% Compute validation predictions
validationPredictions = kfoldPredict(partitionedModel);

% Compute validation RMSE
validationRMSE = sqrt(kfoldLoss(partitionedModel, 'LossFun', 'mse'));
end

function [trainedModel, validationRMSE] = trainEnsembleRegressionPCA(trainingData, responseData)
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData,
% responseData)
% Returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: A matrix with the same number of columns and data type
%       as the matrix imported into the app.
%
%      responseData: A vector with the same data type as the vector
%       imported into the app. The length of responseData and the number of
%       rows of trainingData must be equal.
%
%  Output:
%      trainedModel: A struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: A function to make predictions on new data.
%
%      validationRMSE: A double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input arguments trainingData and responseData.
%
% For example, to retrain a regression model trained with the original data
% set T and response Y, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T, Y)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

% Auto-generated by MATLAB on 10-Nov-2020 17:18:48


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
    'MinLeafSize', 18, ...
    'NumVariablesToSample', 3);
regressionEnsemble = fitrensemble(...
    predictors, ...
    response, ...
    'Method', 'Bag', ...
    'NumLearningCycles', 10, ...
    'Learners', template);

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
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 9 columns because this model was trained using 9 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'});

predictorNames = {'column_1', 'column_2', 'column_3', 'column_4', 'column_5', 'column_6', 'column_7', 'column_8', 'column_9'};
predictors = inputTable(:, predictorNames);
response = responseData(:);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false];

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
        'MinLeafSize', 18, ...
        'NumVariablesToSample', 3);
    regressionEnsemble = fitrensemble(...
        trainingPredictors, ...
        trainingResponse, ...
        'Method', 'Bag', ...
        'NumLearningCycles', 10, ...
        'Learners', template);
    
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


