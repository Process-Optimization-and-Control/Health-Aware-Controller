function [u_,sArray,ehatArray,inputArray,ofValue,solFlag] = SolvingNMPC(solver,x_next,z_next,u_k,mdot_p,GOR,PI,T,nmpcPar)

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
    
    % creating current slack variables
    w0 = [w0;0*ones(nmpcPar.nx,1)];
    lbw = [lbw;0*ones(nmpcPar.nx,1)];
    ubw = [ubw;10000*ones(nmpcPar.nx,1)]; % 10
    
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
    
    % Constraint on erosion
    lbg = [lbg;-inf*ones(nmpcPar.nx,1)]; %zeros
    ubg = [ubg;nmpcPar.x_threshold*ones(nmpcPar.nx,1)];
    
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
    % 10-12: s1
    % 13-15 | 78-80 | 143 - 145: x11, x12, x13
    % 16-77 | 81-142| 146 - 207: z11, z12, z13
    % 208-210: xprev1
    
    % 211-213: u2
    % 214-216: s2
    % 217-219 | 282-284 | 347 - 349: x21, x22, x23
    % 220-281 | 285-346 | 350 - 411: z21, z22, z23
    % 412-414: xprev2
    
    % 415-417: u3
    % 418-420: s3
    % 421-423 | 486-488 | 551 - 553: x31, x32, x33
    % 424-485 | 489-550 | 554 - 615: z31, z32, z33
    % 616-618: xprev3
    
    % total variables in one iteration = 3(u) + 3(s) + 9(xkd) +
    % 186(zkd) + 3(xprev) = 204
    
    u_ = w_opt(7:9);
    
    ofValue = full(sol.f);
    
    ehatArray = [w_opt(1:3), w_opt(208:210)];
    inputArray = w_opt(7:9);
    sArray = w_opt(10:12); %s1
    
    for ii = 1:nmpcPar.np - 1
        temp = 412 + (ii - 1)*204;
        ehatArray = [ehatArray, w_opt(temp:temp + 2)];
        
        temp2 = 214 + (ii - 1)*204;
        sArray = [sArray, w_opt(temp2:temp2 + 2)];
        
        temp3 = 211 + (ii - 1)*204;
        inputArray = [inputArray, w_opt(temp3:temp3 + 2)];
    end

end

end

