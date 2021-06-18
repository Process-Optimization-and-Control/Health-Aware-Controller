clear 
%close all
clc 

%% noise --> For reproducibility
% This is choosing a seed for generating random numbers 
rng(1)  

%% 

%% Initializing the table to store multiple time series:

nTimeseries = 500;
dataSz = [nTimeseries 10];
varNames = {'Num','uArray', 'yMeas','erosionArray', 'H', 'x0', 'xk', 'yk', 'z0', 'zk'};
varTypes = {'double', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
Data = table('Size', dataSz, 'VariableNames', varNames, 'VariableTypes', varTypes);

%% Model initialization
for i_timeseries = 1:nTimeseries

%% Simulation parameters
% Number of simulation steps
nSim = 500;  % time_total = 3600*24*500; %[s]
%sampling time /control interval /1 simulation iteration time
par.T = 3600*24; % [s]

% initial condition (pre-computed)
[x0,z0,u0] = InitialConditionGasLift_5;

% generating sand profile
% 'exp' | 'log' | 'cst' 
% for cst use sand_0 = 0.1
% for log use sand_0 = 0.01, rate = 0.02
% for exp use sand_0 = 0.01, rate = 0.005
%sandArray = sandproductionrate(0.01,500,'log',0.02);
sandArray = sandproductionrate(0.1,500,'cst',0.02);

%% model parameters

par = ParametersGasLift(1,sandArray); % For varying sand production rate 

%states to measurement mapping function
H = zeros(16,length(z0));
%pai - annulus pressure, well 1-3
H(1,1) = 1; 
H(2,2) = 1; 
H(3,3) = 1; 
%pwh - well head pressure, well+ 1-3
H(4,4) = 1; 
H(5,5) = 1; 
H(6,6) = 1; 
%wro - wellhead gas production rate, well 1-3
H(7,25) = 1; 
H(8,26) = 1; 
H(9,27) = 1; 
%wrg - wellhead oil production rate, well 1-3
H(10,28) = 1; 
H(11,29) = 1; 
H(12,30) = 1; 
%prh - riser head pressure
H(13,37) = 1;  
%pm - manifold pressure
H(14,39) = 1;
%wto - riser head total oil production rate
H(15,41) = 1;
%wtg - riser head total gas production rate
H(16,42) = 1; 

%% Simulation initialization
xk = x0;
zk = z0;
uk = u0;
yk = H*z0;

%%
fprintf('Time series number: >>> %0.0f \n',i_timeseries)
%creating random array for the inputs
uArray = [];
%bounds on the inputs
uMin = 1.5;
uMax = 2.5;
for kk = 0:nSim
    if rem(kk,50) == 0 %every 50 days we change the inputs
        uk = (uMax - uMin).*rand(3,1) + uMin;
    end
    uArray = [uArray, uk];
end

% measurements (for plotting)
yMeas = yk;
erosionArray = xk;
for kk = 1:nSim
   
    fprintf('     iteration >>> %0.0f \n',kk)
    
    % integrating the system
    [xk,zk] = WellPlantModel(xk,zk,uArray(:,kk),par);
    par = ParametersGasLift(kk,sandArray);  % Varying sand production 
    par.T = 3600*24;
    
    % Adding noise to the measurements
    yMeas = [yMeas, H*zk + par.scale.*randn(length(yk),1)]; 
    erosionArray = [erosionArray, xk];
end

Data(i_timeseries, :) = {i_timeseries, uArray, yMeas, erosionArray, H, x0, xk, yk, z0, zk};

end

%% saving the data in a mat file 
%filename = 'datamatrix_log_'+string(nTimeseries)+'.mat';
filename = 'datamatrix_cst_'+string(nTimeseries)+'.mat';
%filename = 'datamatrix_'+string(nTimeseries)+'.mat';

save(filename, 'Data')

%% Plotting
figure(1)

time = 0:1:500; %[days]

% System inputs
subplot(2,1,1)
    stairs(time,uArray(1,:),'LineWidth',2);
    hold on
    stairs(time,uArray(2,:),'LineWidth',2);
    stairs(time,uArray(3,:),'LineWidth',2);
    legend('Well 1','Well 2','Well 3');
    
    ylim([0,5]);
    xlabel('Time [day]');
    ylabel('Gas lift rate [kg/s]');

% erosion
subplot(2,1,2);
    plot(time,transpose(Data.erosionArray{1,1}(1,:)),'LineWidth',2);
    hold on
    plot(time,transpose(Data.erosionArray{1,1}(2,:)),'LineWidth',2);
    plot(time,transpose(Data.erosionArray{1,1}(3,:)),'LineWidth',2);
    %ylim([0,5.2]);
    legend('Well 1','Well 2','Well 3');
    legend('Location','northwest');
    xlabel('Time [day]');
    ylabel('Erosion [mm]');

figure(2)

% Pressure
subplot(2,1,1)
    plot(time,yMeas(4,:),'LineWidth',1.5);
    hold on
    plot(time,yMeas(5,:),'LineWidth',1.5);
    plot(time,yMeas(6,:),'LineWidth',1.5);
    axis([0 500 45 60]);
    legend('Well 1','Well 2','Well 3');
    
    %ylim([0,2.2]);
    xlabel('Time [day]');
    ylabel('Well head pressure [bar]');

% erosion
subplot(2,1,2);
    plot(time,yMeas(15,:),'LineWidth',2); %oil
    %plot(time,yMeas(14,:)); %gas
    axis([0 500 50 100]);
    
    xlabel('Time [day]');
    ylabel('Flowrate [kg/s]');
    
    
figure(3)
    % Sand rate profile
    plot(time,sandArray,'o','LineWidth',1.5);
    xlabel('Time [day]');
    ylabel('sand rate [kg/s]');
    title('Sand rate profile (same for all wells)')
