function solver = BuildingNMPC(diff,alg,x_var,z_var,p_var,par,nmpcPar)

import casadi.*

%% Using 3 collocation points:
% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = [0 collocation_points(d, 'legendre')];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:d+1
    if r ~= j
      coeff = conv(coeff, [1, -tau_root(r)]);
      coeff = coeff / (tau_root(j)-tau_root(r));
    end
  end
  % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D(j) = polyval(coeff, 1.0);

  % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);
  for r=1:d+1
    C(j,r) = polyval(pder, tau_root(r));
  end

  % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);
  B(j) = polyval(pint, 1.0);
end

%% Control discretization
h = par.T;

%% Defining system OF
% objective function
w_po = z_var(28:30);

% constraint violation
constV = (x_var - nmpcPar.x_healthy)./(nmpcPar.x_threshold - nmpcPar.x_healthy);

% computing using utility function
L =  -((1 - exp(-constV{1}*w_po{1}))/(constV{1}+0.001) + (1 - exp(-constV{2}*w_po{2}))/(constV{2}+0.001) + (1 - exp(-constV{3}*w_po{3}))/(constV{3}+0.001));% + 1/2 * ((U1 - U_1)'*nmpcPar.R*(U1 - U_1));

% creating system function (LHS of the dynamic equations)
f = Function('f',{x_var,z_var,p_var},{diff/p_var{11},alg,L}); %diff in cm/day --> transforming in cm/s 

%% Defining empty nlp-problem
% objective function
J = 0;

% declare variables (bounds and initial guess)
w = {};
% w0 = [];
% lbw =[];
% ubw = [];

% declare constraints and its bounds
g = {};
% lbg = [];
% ubg = [];

%% declaring parameters
xk_meas = MX.sym('xk_meas',nmpcPar.nx);
zk_meas = MX.sym('zk_meas',nmpcPar.nz);
uk_meas = MX.sym('uk_meas',nmpcPar.nu);

% mdot_p,GOR,PI,T - which are fixed
p = MX.sym('p',8);

%% Lifting initial conditions

% initial state
Xk = MX.sym('X0',nmpcPar.nx);
w = {w{:},Xk};

% initial input
Uk = MX.sym('uk_init',nmpcPar.nu);
w = {w{:}, Uk};

%% Looping through until timeend
for k = 1:nmpcPar.np
    
    % storing the previous input
    uprev = Uk;
    
    % creating current input
    Uk = MX.sym(['U_' num2str(k)],nmpcPar.nu);
    w = {w{:}, Uk}; % 7-9
    
    % Adding constraint for delta_u
    duk = Uk - uprev;
    g = {g{:},duk};
    
    % State at collocation points
    Xkj = {};
    Zkj = {};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)],nmpcPar.nx);
        Zkj{j} = MX.sym(['Z_' num2str(k),'_',num2str(j)],nmpcPar.nz);
        w = {w{:}, Xkj{j},Zkj{j}};
    end
    
    % Loop over collocation points
    Xk_end = D(1)*Xk;
    for j=1:d
        % Expression for the state derivative at the collocation point
        xp = C(1,j+1)*Xk;
        for r=1:d
            xp = xp + C(r+1,j+1)*Xkj{r};
        end
        
        % Append collocation equations
        [fj,gj,qj] = f(Xkj{j},Zkj{j},vertcat(Uk,p));
        g = {g{:}, h*fj - xp};
        g = {g{:}, gj};
        
        % Add contribution to the end state
        Xk_end = Xk_end + D(j+1)*Xkj{j};
        
        % updating objective function
        % Add contribution to quadrature function
        J = J + B(j+1)*qj*h;
    end
    
    % New NLP variable for state at end
    Xk  = MX.sym(['X_',num2str(k+1)],nmpcPar.nx);
    w = {w{:}, Xk};
    
    % Gap
    g = {g{:},Xk_end-Xk};
    
end

% Formalizing problem 
nlp = struct('x',vertcat(w{:}),'g',vertcat(g{:}),'f',J,'p',vertcat(xk_meas,zk_meas,uk_meas,p));

% Assigning solver (IPOPT)
solver = nlpsol('solver','ipopt',nlp);


end

