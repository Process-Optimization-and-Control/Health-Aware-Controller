function [u_,ehatArray,inputArray,ofValue,solFlag] = SolvingNMPC(solver,x_next,z_next,u_k,mdot_p,GOR,PI,T,nmpcPar)

% declare variables (bounds and initial guess)
w0 = [];
lbw =[];
ubw = [];

% declare constraints and its bounds
lbg = [];
ubg = [];

% initial state
lbw = [lbw,x_next]; 
ubw = [ubw,x_next];
w0 = [w0;x_next];

% initial input
w0 = [w0;u_k];
lbw = [lbw;u_k];
ubw = [ubw;u_k];

%% Looping through until timeend
for k = 1:nmpcPar.np
    w0 = [w0; u_k];
    lbw = [lbw;nmpcPar.umin*ones(nmpcPar.nu,1)];
    ubw = [ubw;nmpcPar.umax*ones(nmpcPar.nu,1)];
    
%     % creating current slack variables
%     w0 = [w0;0*ones(nmpcPar.nx,1)];
%     lbw = [lbw;0*ones(nmpcPar.nx,1)];
%     ubw = [ubw;10000*ones(nmpcPar.nx,1)]; % 10
    
    if k > nmpcPar.nm
        lbg = [lbg;zeros(nmpcPar.nu,1)];
        ubg = [ubg; zeros(nmpcPar.nu,1)];
    else
        lbg = [lbg;-nmpcPar.dumax*ones(nmpcPar.nu,1)];
        ubg = [ubg;nmpcPar.dumax*ones(nmpcPar.nu,1)];
    end
    
    for d = 1:3
        % creating states at collocation points
        w0 = [w0;x_next;z_next];
        lbw = [lbw;zeros(nmpcPar.nx,1);zeros(nmpcPar.nz,1)];
        ubw = [ubw;10*ones(nmpcPar.nx,1);inf*ones(nmpcPar.nz,1)];%inf
     
    end
    
    % integrating the system
    for d = 1:3 
     
        % Adding xk and Xk1 as constrains as they must be equal - in
        % collocation intervals
        % algebraic constraints are set to zero in the collocation point
        lbg = [lbg;zeros(nmpcPar.nx,1);zeros(nmpcPar.nz,1)];
        ubg = [ubg;zeros(nmpcPar.nx,1);zeros(nmpcPar.nz,1)];  
    end
    
       
    % New NLP variable for state at end
    w0 = [w0;x_next];
    lbw = [lbw;zeros(nmpcPar.nx,1)];
    ubw = [ubw;10*ones(nmpcPar.nx,1)]; %inf
    
    % Gap
    lbg = [lbg;zeros(nmpcPar.nx,1)];
    ubg = [ubg;zeros(nmpcPar.nx,1)];
    
%     % Constraint on erosion
%     lbg = [lbg;-inf*ones(nmpcPar.nx,1)]; %zeros
%     ubg = [ubg;nmpcPar.x_threshold*ones(nmpcPar.nx,1)];
    
end


% Solving the problem
sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg,'p',[x_next;z_next;u_k;mdot_p;GOR;PI;T]);

%% Extracting solution
w_opt = full(sol.x);

% catch error
if solver.stats.success ~=1
    % solution failed 
    solFlag = 0;
    
    u_ = u_k; % fix inputs

    %dummy
    sArray = zeros(nmpcConfig.nx,1);
    ehatArray = zeros(nmpcConfig.nx,nmpcConfig.np + 1);

else
    % solution succeeded 
    solFlag = 1;
    
    % variable order  
    % 1-3: x0
    % 4-6: u0
    % 7-9: u1
    % 10-12 | 75-77 | 140 - 142: x11, x12, x13
    % 13-74 | 78-139| 143 - 204: z11, z12, z13
    % 205-207: xprev1
    
    % 208-210: u2
    % 211-213 | 276-278 | 341 - 343: x21, x22, x23
    % 214-275 | 279-340 | 344 - 405: z21, z22, z23
    % 406-408: xprev2
    
    % 409-411: u3
    % 412-414 | 477-479 | 542 - 544: x31, x32, x33
    % 415-476 | 480-541 | 545 - 606: z31, z32, z33
    % 607-609: xprev3
    
    % total variables in one iteration = 3(u) + 9(xkd) +
    % 186(zkd) + 3(xprev) = 201
    
    u_ = w_opt(7:9);
    
    ofValue = full(sol.f);
    
    ehatArray = [w_opt(1:3), w_opt(205:207)];
    inputArray = w_opt(7:9);
    
    for ii = 1:nmpcPar.np - 1
        temp = 406 + (ii - 1)*201;
        ehatArray = [ehatArray, w_opt(temp:temp + 2)];
               
        temp3 = 208 + (ii - 1)*201;
        inputArray = [inputArray, w_opt(temp3:temp3 + 2)];
    end

end

end

