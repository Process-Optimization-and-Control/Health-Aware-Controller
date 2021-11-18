% Implementation of Health Aware Controller based on Joachim's thesis
clear 
close all
clc 

import casadi.*

% for reproducibility
rng('default');

%% simulation configuration
nMC = 200; 
simLength = 300; % 300 | 450 | 500 | 550

% Initial condition
[x0,z0,u0] = InitialConditionGasLift_5;

% System parameters
par = ParametersGasLift;

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
par.H = H;

% generating sand profile
sandArray = sandproductionrate(0.1,simLength,'cst',0.02);

% tuning (from HAC)
nmpcConfig.dumax = 0.05; % 0.01 | 0.1
nmpcConfig.umax = 2;
nmpcConfig.umin = 0.4;

%% generating input sequence and initial erosion 
% sampling the inputs
eps_MC = 0.01*exprnd(1, [nMC, 3]);
%histogram(eps_MC(:,1),20)

% sampling the initial time of the controller actions
t_MC = 150*wblrnd(1,10,[nMC, 3]);
%histogram(t_MC(:,1),20)

% sampling the delta u's
for well = 1:3
    deltaU_MC{well} = 0 + (0 - nmpcConfig.dumax)*betarnd(1,5,[nMC,simLength]);
    %histogram(deltaU_MC{well},20)
end

% building input sequence (always starts with nmpcConfig.umax)
for well = 1:3
    U_MC{well} = [nmpcConfig.umax*ones(nMC,1), zeros(nMC,simLength - 1)];
end

for ii = 1:nMC
    % "running the simulation
    for kk = 2:simLength
        temp = [];
        for well = 1:3
            % checking if the MPC action time has started
            if kk <= t_MC(ii,well)
                temp1 = U_MC{well}(ii,kk - 1);
            else
                % if yes, add DU
                temp1 = U_MC{well}(ii,kk - 1) + deltaU_MC{well}(ii,kk);
            
                % enforce bounds
                if temp1 < nmpcConfig.umin
                    temp1 = nmpcConfig.umin;
                end
            end
            
            % update input
            U_MC{well}(ii,kk) = temp1;
        end
    end
end

% % plotting inputs
% figure(1)
% for well = 1:3
%     subplot(3,1,well)
%     hold on
%     for ii = 1:nMC
%        plot(1:simLength,U_MC{well}(ii,:),'Color',[0,0,0,0.05],'Linewidth',1)
%     end
%     
%     ylim([nmpcConfig.umin nmpcConfig.umax])
%     ylabel('Q_g [L/min]','FontSize',10)
%     
%     xticks(0:50:(simLength - 1))
%     xlim([0 (simLength - 1)])
%     xlabel('time [days]','FontSize',10)
%     
%  
% end

%% running the model
% m_k
% phenomeno = 0;
% stepwise = 0;
% neuralnet = 1;
flagModel = zeros(3,1); %[phenomeno,stepwise,neuralnet];

for m_k = 1:3
    for well = 1:3
        erosionArray{m_k,well} = [];
    end
end

% Monte carlo simulation
for m_k = 1:3
    flagModel = zeros(3,1); %[phenomeno,stepwise,neuralnet];
    flagModel(m_k) = 1;
    
    for ii = 1:nMC
       %% Initializing the simulation
        xk = eps_MC(ii,:)'; % from the sampled initial condition
        zk = z0;
        uk = [U_MC{1}(ii,1); U_MC{2}(ii,1); U_MC{3}(ii,1)];
    
        xPlant = xk;
        for kk = 2:simLength
            fprintf(' model: %0.0f | MC iteration: %0.0f | simulation time: %0.0f [day]\n',m_k,ii,kk)
            % picking up from the 
            uTemp = [U_MC{1}(ii,kk); U_MC{2}(ii,kk); U_MC{3}(ii,kk)];
            
            [xk,zk] = WellPlantModel(xk,zk,uTemp,sandArray(kk + 1),par,flagModel);
            xPlant = [xPlant, xk];
            
        end
        for well = 1:3
            erosionArray{m_k,well} = [erosionArray{m_k,well}; xPlant(well,:)];
        end
    end
end

name = 'openLoopAnalysis';
save(name,'simLength','erosionArray','U_MC','deltaU_MC','eps_MC','t_MC')

%% plotting degradation
figure(2)
for m_k = 1:3
    subplot(3,1,m_k)
    hold on
    for ii = 1:nMC
       plot(1:simLength,erosionArray{m_k,1},'Color',[0,0,0,0.05],'Linewidth',1)
    end
    
    %ylim([0 nmpcConfig.umax])
    ylabel('Crack length [cm]','FontSize',10)
    
    xticks(0:50:(simLength - 1))
    xlim([0 (simLength - 1)])
    xlabel('time [days]','FontSize',10)
    

end

