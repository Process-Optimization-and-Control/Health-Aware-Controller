% Implementation of Health Aware Controller based on Joachim's thesis
clear 
close all
clc 

import casadi.*

% for reproduciability
rng('default');

%% Configuration
% Flag to choose the model (choose only one)
phenomeno = 1;
stepwise = 0;
neuralnet = 0;
flagModel = [phenomeno,stepwise,neuralnet];

% show the plot while simulating
showPlot = false; % true | false

% setting simulation length
simLength = 500; % 300 | 450 | 500 | 550

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

% Building Dynamic Model
[diff,alg,x_var,z_var,p_var] = BuildingDynModel(par,flagModel);

% control tuning
nmpcConfig.umax = 2;
nmpcConfig.umin = 0.4;
nmpcConfig.dumax = 0.05; % 0.01 | 0.1

%crack length region
nmpcConfig.x_healthy = 0.1; 
nmpcConfig.x_threshold = 0.65;   % 0.4 | 0.65
nmpcConfig.x_broken = 0.6;      % 0.5 | 0.75
nmpcConfig.nx = size(x0,1);
nmpcConfig.nz = size(z0,1);
nmpcConfig.nu = size(u0,1);

nmpcConfig.nm = 20;  % 70
nmpcConfig.np = 50; % 100
nmpcConfig.rho = 1e3*eye(nmpcConfig.nx); % 1e3 | 999999
nmpcConfig.R = 1*eye(nmpcConfig.nu); % 0 | 0.01 | 1  

% Building NMPC
solverNMPC = BuildingNMPC(diff,alg,x_var,z_var,p_var,par,nmpcConfig);

%% Initializing the simulation
xk = x0;
xHat = x0; % initial value is assumed known 
zk = z0;
uk = u0;
yk = par.H*z0;

%initializing slacks (not needed - just for plotting if not using the controller)
s = zeros(nmpcConfig.nx,nmpcConfig.np); 
erosionHat = zeros(nmpcConfig.nx,nmpcConfig.np + 1); % erosion predict "inside the controller" --> prognostics
inputSeq = zeros(nmpcConfig.nx,nmpcConfig.np); % input sequence prediction

OF = 0;

solFlag = 0;
controlTime = 0; 

% creating variable to save breakdown time (if applicable)
tBreak = [];

% preparing for plotting
%colors associated with each well
cc = {'b','k','r'};

% plant info
xPlant = [];
zPlant = [];
yPlant = [];
totalProductionPlant = [];

%diagnositics info
erosionHatDiag = [];

%control info
inputControl = [];
flagControl = [];
ofControl = [];
CPUtimeControl = [];

for ii = 0:simLength
    erosionPrediction{ii + 1} = [];
    inputPrediction{ii + 1} = [];
end

%% Simulating

for tt = 0:simLength
    fprintf('     iteration >>> %0.0f [day]\n',tt)

    if tt ~=0
        % estimating the current erosion
        xHat = Diagnostics(xk,xHat,yk,uk,sandArray(tt),par,flagModel);
        
        tic
        % Calculating input with NMPC
        % perfect information
        [uk,erosionHat,inputSeq,OF,solFlag] = SolvingNMPC(solverNMPC,xk,zk,uk,sandArray(tt),par.GOR,par.PI,par.T,nmpcConfig);

        % with diagnostics
        % [uk,s,erosionHat,inputSeq,OF,solFlag] = SolvingNMPC(solverNMPC,xHat,zk,uk,sandArray(tt),par.GOR,par.PI,par.T,nmpcConfig);

          % computing execution time
          controlTime = toc;
    end
    
    % Finding the state after applying input u
    [xk,zk] = WellPlantModel(xk,zk,uk,sandArray(tt + 1),par);
    
        %%%%%%%%%%%%%%%
        % Saving Data %
        %%%%%%%%%%%%%%%
        erosionHatDiag = [erosionHatDiag, xHat];

        inputControl = [inputControl, uk];
        flagControl =  [flagControl, solFlag];
        erosionPrediction{tt + 1} = erosionHat;
        inputPrediction{tt + 1} = inputSeq;
        ofControl = [ofControl, OF];
        CPUtimeControl = [CPUtimeControl, controlTime];

        xPlant = [xPlant, xk];
        zPlant = [zPlant, zk];

        yk = par.H*zk + par.scale.*randn(length(yk),1);
        yPlant = [yPlant, yk]; 
        totalProductionPlant = [totalProductionPlant, sum(zk(28:30))];
    
      

        if tt > 0 && showPlot %rem(t,50) == 0 %checking results sporadiacally
            time = 0:1:tt;
             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plotting --- controller %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            f1 = figure(1);
            clf
            
            subplot(3,1,1);
                hold on
                for well = 1:3
                    stairs(time,inputControl(well,:),cc{well},'LineWidth',1.5,'HandleVisibility','off');
                end 

                xline(simLength,'r:','HandleVisibility','off');
                yline(nmpcConfig.umax,'k--','LineWidth',1);

                ylim([0,2.75]);
                xlim([0,simLength + 100]);

                legend('Max. gas cap.','Position',[0.69 0.90 0.20 0.041])

                xlabel('time [day]');
                ylabel('Gas lift rate [kg/s]');
            
            subplot(3,1,2);

                yline(nmpcConfig.x_threshold,'k--','LineWidth',1);
                hold on
                for well = 1:3
                    plot(time,xPlant(well,:),cc{well},'LineWidth',1.5);
                    plot(tt:tt + nmpcConfig.np,erosionHat(well,:),cc{well},'linestyle',':','LineWidth',1);
                    plot(tt,xHat(well),'MarkerFaceColor',cc{well},'marker','o','HandleVisibility','off');
                end 

                xline(simLength,'r:');

                text(simLength + 50,0.7,'\leftarrow Maintenance','HorizontalAlignment','center','FontSize',7);
                text(simLength + 50,0.5,'Stop','HorizontalAlignment','center','FontSize',7);

                legend({'Limit','Real','Predicted'},'Position',[0.14 0.59 0.17 0.11])

                box on

                ylim([0,nmpcConfig.x_threshold*1.5]);
                xlim([0,simLength + 100]);
                xlabel('Time [day]');
                ylabel('Erosion [mm]');

            subplot(3,1,3);
                plot(time,totalProductionPlant,'b','LineWidth',1.5);
                hold on
                xline(simLength,'r:');

                ylim([80,88]);
                xlim([0,simLength + 100]);
                xlabel('Time [day]');
                ylabel('Oil production [kg/s]');
                     
        end
        
        % emulating system breaking if threshold is achieved
        if ~all(xk < nmpcConfig.x_broken)
           
            % saving break time
            tBreak = tt;
            
            % breaking loop
            break 
            
            beep
        end
        
    
end

%% saving results
if flagModel(1) == 1
    name = 'HAC_cst_pheno';
elseif flagModel(2) == 1
    name = 'HAC_cst_step';  
elseif flagModel(3) == 1
    name = 'HAC_cst_NN';
end

save(name,'simLength','tBreak','xPlant','zPlant','yPlant','totalProductionPlant','inputControl','flagControl','erosionPrediction','CPUtimeControl','ofControl','erosionHatDiag','inputPrediction')

%% Plotting the final profiles
f2 = figure(2);

if ~isempty(tBreak)
   simLength = tBreak;
end
time = 0:1:simLength;

    subplot(3,1,1);
        stairs(time,inputControl(1,:),cc{1},'LineWidth',1.5,'HandleVisibility','off');
        hold on 
        stairs(time,inputControl(2,:),cc{2},'LineWidth',1.5,'HandleVisibility','off');
        stairs(time,inputControl(3,:),cc{3},'LineWidth',1.5,'HandleVisibility','off');
        
        yline(nmpcConfig.umax,'k--','LineWidth',1);
        
        ylim([0,2.75]);
        xlim([0,simLength]);

        legend('Max. gas cap.','Position',[0.69 0.90 0.20 0.041])
    
        xlabel('time [day]');
        ylabel('Gas lift rate [kg/s]');

    subplot(3,1,2);

        yline(nmpcConfig.x_threshold,'k--','LineWidth',1);
        hold on
        plot(time,xPlant(1,:),cc{1},'LineWidth',1.5);
        plot(time,xPlant(2,:),cc{2},'LineWidth',1.5);
        plot(time,xPlant(3,:),cc{3},'LineWidth',1.5);
                   
        box on
        
        legend({'limit','Well 1','Well 2','Well 3'},'Location','northwest');
        
        ylim([0,3]);
        xlim([0,simLength]);        
        xlabel('Time [day]');
        ylabel('Erosion [mm]');

    subplot(3,1,3);
        plot(time,totalProductionPlant,'b','LineWidth',1.5);
        hold on
        
        ylim([80,88])
        xlim([0,simLength]);
        xlabel('Time [day]');
        ylabel('Total oil production [kg/s]');
    
    set(gcf,'color','w');

    namePlot = [name,'_results.pdf'];
    print(f2,'-r300','-dpdf',namePlot);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting --- diagnosis  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    f3 = figure(3);    
    
    for well = 1:3
        subplot(3,1,well);
        yline(nmpcConfig.x_threshold,'k--','LineWidth',1);
        hold on
        plot(time,erosionHatDiag(well,:),cc{well},'LineWidth',1.5);
        plot(time,xPlant(well,:),cc{well},'linestyle',':','LineWidth',1.5);
        
        box on
        
        legend({'Limit','Estimated','Real'},'Location','northwest');
        
        ylim([0,nmpcConfig.x_threshold*1.5]);
        xlim([0,simLength]);
        xlabel('Time [day]');
        ylabel('Erosion [mm]');
        
        set(gcf,'color','w');
        
    end

    namePlot = [name,'_diag_results.pdf'];
    print(f3,'-r300','-dpdf',namePlot);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting --- statistics %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    f4 = figure(4);    
    
  subplot(2,1,1)
        plot(time,flagControl,'marker','x','linestyle',':','markersize',5);
        
        ylim([-0.1,1.1]);
        yticks(0:1);
        yticklabels({'no','yes'});
        xlim([0,simLength]);
        xlabel('iteration [-]');
        title('Control converged?');
       
       
  subplot(2,1,2)
        plot(time,CPUtimeControl,'marker','x','linestyle',':','markersize',5);
        
        xlim([0,simLength]);
        xlabel('Iteration [-]');
        ylabel('exec. time [s]');
        title('Control time');
  
set(gcf,'color','w');

    namePlot = [name,'_stats.pdf'];
    print(f4,'-r300','-dpdf',namePlot);
 
%     %%%%%%%%%%%%%%%%%%%%%%%
%     % Plotting --- slacks %
%     %%%%%%%%%%%%%%%%%%%%%%%
%     f5 = figure(5);    
%     
%     cmap = flip(gray(simLength));
%     
%     for well = 1:3
%         subplot(3,1,well)
%         hold on
%         for ii = 1:simLength
%             plot(ii:(ii + nmpcConfig.np - 1),slackControl{1,ii}(well,:),'Color',cmap(ii,:),'LineWidth',0.75);
%         end
%         
%         box on
%         
%         ylim([0,1]);
%         xlim([0,simLength + nmpcConfig.np]);
%         xticks(0:50:500);
%         xlabel('Iteration [-]');
%         ylabel('slack [mm]');
%         title(['Well ',num2str(well)]);
%         
%     end
% 
%                   
% set(gcf,'color','w');
% 
%     namePlot = [name,'_slack.pdf'];
%     print(f5,'-bestfit','-r300','-dpdf',namePlot);
%     
    %%%%%%%%%%%%%%%%%%%%%%%
    % Plotting --- inputs %
    %%%%%%%%%%%%%%%%%%%%%%%
    f6 = figure(6);    
    
    cmap = flip(gray(simLength));
    
    for well = 1:3
        subplot(3,1,well)
        hold on
        for ii = 1:simLength
            plot(ii:(ii + nmpcConfig.np - 1),inputPrediction{1,ii}(well,:),'Color',cmap(ii,:),'LineWidth',0.75);
        end
        
        yline(nmpcConfig.umax,'k--','LineWidth',1);

        box on
        
        ylim([0,2.75]);
        xlim([0,simLength]);
    
        xlabel('Iteration [-]');
        ylabel('Gas lift rate [kg/s]');
        title(['Well ',num2str(well)]);
        
    end

                  
set(gcf,'color','w');

    namePlot = [name,'_input_seq.pdf'];
    print(f6,'-bestfit','-r300','-dpdf',namePlot);

%    %% 
%     %%%%%%%%%%%%%%%%%%%%%%%
%     % Plotting --- OF     %
%     %%%%%%%%%%%%%%%%%%%%%%%
%     f7 = figure(7);    
%     
%     % computing terms
%     RtermArray = [];
%     StermArray = [];
%     for ii = 1:simLength
%     % Slack
%         temp1 = slackControl{1,ii};
%     % Regularization
%         temp2 = inputPrediction{1,ii} - [u0, inputPrediction{1,ii}(:,1:end - 1)];
%         
%         Rterm = 0;
%         Sterm = 0;
%         for kk = 1:nmpcConfig.np
%             Rterm = Rterm + temp2(:,kk)'*nmpcConfig.R*temp2(:,kk);
%             Sterm = Sterm + nmpcConfig.rho*sum(temp1(:,kk));
%         end
%             
%          RtermArray = [RtermArray, Rterm];
%          StermArray = [StermArray, Sterm];
%     end
%     
%         plot(time,ofControl,'k','marker','d','linestyle','-','markersize',3);
%         hold on 
%         plot(time(2:end),RtermArray,'r','marker','x','linestyle',':','markersize',3);
%         plot(time(2:end),StermArray,'b','marker','o','linestyle','--','markersize',3);
% 
%         legend({'OF value','Reg. Term','Slack Term'})
%         
%         xlim([0,simLength]);
%         xlabel('iteration [-]');
%         ylabel('OF value [-]');
%     
%                   
% set(gcf,'color','w');
% 
%     namePlot = [name,'_obj_fun.pdf'];
%     print(f7,'-bestfit','-r300','-dpdf',namePlot);  
%     

