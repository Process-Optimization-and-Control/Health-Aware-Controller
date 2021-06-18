function [diff,alg,x_var,z_var,p_var] = BuildingDynModel(par,modelFlag)

    import casadi.*

    %% normalization (for stepwise and NN)
    norm = load('normalizationValues');
    
    %% Parameters
    %number of wells
    n_w = 3; %[]
    %gas constant
    R = par.R; %[m3 Pa /K /mol]
    %molecular weigth
    Mw = par.Mw; %[kg/mol]

    %properties
    %density of oil - dim:  nwells x 1
    rho_o = par.rho_o; %[kg/m3]
    %riser oil density 
    rho_ro = par.rho_ro;%[kg/m3]
    %1cP oil viscosity
    mu_oil = par.mu_oil;% [Pa s] 

    %project
    %well parameters - dim:  nwells x 1
    L_w = par.L_w; %[m]
    H_w = par.H_w; %[m]
    D_w = par.D_w; %[m]
    A_w = par.A_w;%[m2]

    %well below injection - [m]
    L_bh = par.L_bh;
    H_bh = par.H_bh;
    D_bh = par.D_bh;
    A_bh = par.A_bh;%[m2]

    %annulus - [m]
    H_a = par.H_a;
    V_a = par.V_a; %[m3]

    %riser - [m]
    L_r = par.L_r;
    H_r = par.H_r;
    D_r = par.D_r;
    A_r = par.A_r;%[m2]

    %injection valve characteristics  - dim:  nwells x 1
    C_iv = par.C_iv;%[m2]
    %production valve characteristics  - dim:  nwells x 1
    C_pc = par.C_pc;%[m2]
    %riser valve characteristics
    C_pr = par.C_pr;%[m2]
    % account for differences in the vapor and oil velocity
    slip = par.slip_real;


    %% For erosion model
    % Sand
    d_p = par.d_p; %[m] particle diameter
    % Choke
    D = par.D; %[m] Gap between body and cage
    % Constants
    n = par.n; %[-] Velocity coefficient
    % Precalculations of erosion in choke:
    A_g = par.A_g; %[m2] Effective gallery area
    gma = d_p./D;

    %% Algebraic states
    %pressure - annulus
    p_ai = MX.sym('p_ai',n_w);      % 1-3 [bar] (bar to Pa  = x10^5)
    %pressure - well head
    p_wh = MX.sym('p_wh',n_w);      % 4-6 [bar]
    %pressure - injection point
    p_wi = MX.sym('p_wi',n_w);      % 7-9 [bar]
    %pressure - below injection point (bottom hole)
    p_bh = MX.sym('p_bh',n_w);      % 10-12 [bar]
    %density - annulus
    rho_ai = MX.sym('rho_ai',n_w);  % 13-15 [100 kg/m3]
    %mixture density in tubing
    rho_m = MX.sym('rho_m',n_w);    % 16-18 [100 kg/m3]
    %well injection flow rate
    w_iv = MX.sym('w_iv',n_w);      % 19-21 [kg/s]
    %wellhead total production rate
    w_pc = MX.sym('w_pc',n_w);      % 22-24 [kg/s]
    %wellhead gas production rate
    w_pg = MX.sym('w_pg',n_w);      % 25-27 [kg/s]
    %wellhead oil production rate
    w_po = MX.sym('w_po',n_w);      % 28-30 [kg/s]
    %oil rate from reservoir
    w_ro = MX.sym('w_ro',n_w);      % 31-33 [kg/s]
    %gas rate from reservoir
    w_rg = MX.sym('w_rg',n_w);      %34 -36 [0.1 kg/s]
    %riser head pressure
    p_rh = MX.sym('p_rh',1);        % 37 [bar]
    %mixture density in riser
    rho_r =  MX.sym('rho_r',1);     % 38 [100 kg/s]
    %manifold pressure
    p_m = MX.sym('p_m',1);          % 39 [bar]
    %riser head total production rate
    w_pr = MX.sym('w_pr',1);        % 40 [kg/s]
    %riser head total oil production rate
    w_to = MX.sym('w_to',1);        % 41 [kg/s]
    %riser head total gas production rate
    w_tg = MX.sym('w_tg',1);        % 42 [kg/s]
    %gas holdup @ annulus 
    m_ga = MX.sym('m_ga',n_w); % 43-45 [ton]
    %gas holdup @ well
    m_gt = MX.sym('m_gt',n_w); % 46-48 [ton]
    %oil holdup @ well 
    m_ot = MX.sym('m_ot',n_w); % 49-51 [ton]
    %gas holdup @ riser
    m_gr = MX.sym('m_gr',1);   % 52 [ton]
    %oil holdup @ riser
    m_or = MX.sym('m_or',1);   % 53 [ton]
    %particle impact velocity
    V_p = MX.sym('V_p',n_w); % 54-56
    %dynamic viscosity of mixture
    mu_f = MX.sym('mu_f',n_w); % 57-59
    g1 = MX.sym('g1',n_w); % 60-62

    %control input
    %gas lift rate
    w_gl = MX.sym('w_gl',n_w); %[kg/s]

    %parameters
    p_res = MX.sym('p_res',n_w);
    %productivity index
    PI = MX.sym('PI',n_w); %[kg s^-1 bar-1]
    %GasOil ratio
    GOR = MX.sym('GOR',n_w); %[kg/kg]
    %Annulus temperature
    T_a = MX.sym('T_a',n_w); %[oC]
    %well temperature
    T_w = MX.sym('T_w',n_w); %[oC]
    %riser temperature
    T_r = MX.sym('T_r',1); %[oC]
    %separator pressure
    p_s = MX.sym('p_s',1); %[bar]
    %time transformation: CASADI integrates always from 0 to 1 and the USER does the time
    %scaling with T.
    T = MX.sym('T',1); %[s]
    
    %sandrate rate
    mdot_p = MX.sym('mdot_p',1); %

    %erosion rate
    ER = MX.sym('ER',n_w); %

    %% Modeling
    %gas fraction (mass) of the well holdup - avoiding zero division
    xGwH = (m_gt.*1e3./max(1e-3,(m_gt.*1e3+m_ot.*1e3)));
    %gas fraction (mass) of the riser holdup
    xGrH = (m_gr.*1e3./(m_gr.*1e3+m_or.*1e3));

    xGw = slip.*xGwH./(1 + (slip-1).*xGwH);
    xOw = 1 - xGw;
    xGr = slip.*xGrH./(1 + (slip-1).*xGrH);
    xOr = 1 - xGr;

    % ===================================
    %     Well model with/withou pressure loss
    % ===================================
    % algebraic equations (all symbolic)
    %annulus pressure - %g = 9.81
    f1 = -p_ai.*1e5 + ((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3) + (Mw./(R.*T_a).*((R.*T_a./(V_a.*Mw) + 9.81.*H_a./V_a).*m_ga.*1e3)).*9.81.*H_a;
    %well head pressure
    f2 = -p_wh.*1e5 + ((R.*T_w./Mw).*(m_gt.*1e3./(L_w.*A_w + L_bh.*A_bh - m_ot.*1e3./rho_o))) - ((m_gt.*1e3+m_ot.*1e3 )./(L_w.*A_w)).*9.81.*H_w/2;
    %well injection point pressure
    f3 = -p_wi.*1e5 + (p_wh.*1e5 + 9.81./(A_w.*L_w).*max(0,(m_ot.*1e3+m_gt.*1e3-rho_o.*L_bh.*A_bh)).*H_w) + (128.*mu_oil.*L_w.*w_pc./(3.14.*D_w.^4.*((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3)));
    %bottom hole pressure
    f4 = -p_bh.*1e5 + (p_wi.*1e5 + rho_o.*9.81.*H_bh + 128.*mu_oil.*L_bh.*w_ro./(3.14.*D_bh.^4.*rho_o));
    %gas density in annulus
    f5 = -rho_ai.*1e2 +(Mw./(R.*T_a).*p_ai.*1e5);
    %fluid mixture density in well
    f6 = -rho_m.*1e2 + ((m_gt.*1e3 + m_ot.*1e3).*p_wh.*1e5.*Mw.*rho_o)./(m_ot.*1e3.*p_wh.*1e5.*Mw + rho_o.*R.*T_w.*m_gt.*1e3);
    %well injection flow rate
    f7 = -w_iv + C_iv.*sqrt(rho_ai.*1e2.*max(0,(p_ai.*1e5 - p_wi.*1e5)));
    %wellhead prodution rate
    f8 = -w_pc + 1.*C_pc.*sqrt(rho_m.*1e2.*max(0,(p_wh.*1e5 - p_m.*1e5)));
    %wellhead gas production rate
    f9 = -w_pg + xGw.*w_pc;
    %wellhead oil prodution rate
    f10 = -w_po + xOw.*w_pc;
    %oil from reservoir flowrate
    f11 = -w_ro + PI.*1e-6.*(p_res.*1e5 - p_bh.*1e5);
    %gas from reservoir production rate
    f12 = -w_rg.*1e-1 + GOR.*w_ro;  
    %riser head pressure
    f13 = -p_rh.*1e5 + ((R.*T_r./Mw).*(m_gr.*1e3./(L_r.*A_r))) - ((m_gr.*1e3+m_or.*1e3 )./(L_r.*A_r)).*9.81.*H_r/2;
    %riser density
    f14 = -rho_r.*1e2 + ((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*T_r.*m_gr.*1e3);
    %manifold pressure
    f15 = -p_m.*1e5 + (p_rh.*1e5 + 9.81./(A_r.*L_r).*(m_or.*1e3+m_gr.*1e3).*H_r) + (128.*mu_oil.*L_r.*w_pr./(3.14.*D_r.^4.*((m_gr.*1e3 + m_or.*1e3).*p_rh.*1e5.*Mw.*rho_ro)./(m_or.*1e3.*p_rh.*1e5.*Mw + rho_ro.*R.*T_r.*m_gr.*1e3)));
    %total production rate of well
    f16 = -w_pr + 1.*C_pr.*sqrt(rho_r.*1e2.*(p_rh.*1e5 - p_s.*1e5));
    %oil total production rate
    f17 = -w_to + xOr.*w_pr; 
    %gas total production rate
    f18 = -w_tg + xGr.*w_pr; 
    % setting differential equations as algebraic equations since the
    % dynamics of ER is on a much larger time scale
    f19 = (w_gl - w_iv).*1e-3;
    f20 = (w_iv + w_rg.*1e-1 - w_pg).*1e-3;
    f21 = (w_ro - w_po).*1e-3;
    f22 = (sum(w_pg) - w_tg).*1e-3 ;
    f23 = (sum(w_po) - w_to).*1e-3 ;
    f24 = - V_p + 3/(4*A_g).*(w_po./rho_o + R*T_w.*w_pg./(p_wh.*10^5*Mw));
    f25 = - mu_f + mu_oil.*(w_po./rho_o)./(w_po./rho_o + R.*T_w.*w_pg./(p_wh.*10^5*Mw));
    f26 = -g1 + gma/0.1;
    % differential equations - (all symbolic) - [ton]
    % Erosion rate

    if modelFlag(1) == 1
        % Use phenomenological model
        ER_constant = par.K*par.F*par.C_1*par.GF*mdot_p*par.C_unit/(par.rho_t*par.A_t);
        df1 =  ER_constant.*g1.*(V_p).^n;
    else
        % computing the regressors
        %regr = [mdot_p; p_ai; p_wh; w_ro; w_rg; p_rh; p_m; w_to; w_tg; w_gl];
        % separating the regressor for the three wells. Some regressors are
        % shared by the three models
        regr = [[mdot_p,mdot_p,mdot_p]; 
                p_ai'; 
                p_wh'; 
                w_ro'; 
                w_rg'; 
                [p_rh,p_rh,p_rh];
                [p_m,p_m,p_m];
                [w_to,w_to,w_to];  
                [w_tg,w_tg,w_tg]; 
                w_gl'];        

        % normalizing regressors with data from (ExpDataAnalysis.m)
        regrN = (regr - norm.regrCenter)./norm.regrScale;
        
        if modelFlag(2) == 1
            % number of regressors
            nReg = length(regr);
            
            % modeling stepwise 
            % putting the weights in the right place
            SW_M1 = zeros(1,nReg);
            SW_M1(1) =  0.936202622379369;
            SW_M1(2) = -0.00210657217345035;
            SW_M1(3) =  0.00340292568131542;
            SW_M1(4) = -0.00201748541907079;
            SW_M1(5) =  0.00805275248024692;
            SW_M1(8) =  0.00991932368893031;
            SW_M1(9) =  0.00896809727051577;
            SW_M1(10) = 0.102758189995450;
            
            SW_M2 = zeros(nReg,nReg);
            SW_M2(1,2) =  0.0324197329838858;
            SW_M2(1,10) = 0.0780912090793351;
            SW_M2(2,4) = -0.00874738623206122;
            SW_M2(2,5) = -0.0258263927628834;
            SW_M2(2,8) =  0.160177833500558;
            SW_M2(2,9) =  0.00893152554186846;
            SW_M2(2,10) = 0.0117711262910581;
            SW_M2(5,8) = -0.126740929538979;
            SW_M2(5,10) = 0.0281549336959051;
            SW_M2(8,9) = -0.0606654291484604;
            
            intercept = -0.00741542407032436;            

            %building model
            df1N = [];
            for well = 1:3
                df1N = [df1N, intercept + SW_M1*regrN(:,well) + regrN(:,well)'*SW_M2*regrN(:,well)];
            end
            
        elseif modelFlag(3) == 1
            % Use neural net model
            % Input 1
            x1_step1.xoffset = [-0.960381876486174;-4.46248802722717;-5.17852887458485;-4.48073715332605;-3.39270547151972;-4.89034267011056;-4.09638956423933;-22.1011864634595;-12.6725101583375;-1.68741224378933];
            x1_step1.gain = [0.744587117901481;0.305647520874011;0.217139586607698;0.278754220747781;0.313094454989053;0.224952360464966;0.2398799497799;0.0883066374206199;0.130806471865733;0.582950679298246];
            x1_step1.ymin = -1;

            % Layer 1
            b1 = [0.14367720318601545637;2.6585932597481978235;0.21614754193338026056;-0.77521005365754414029;-0.15016048001209872376;-0.18895893794020338086;0.33334021415001674482;1.0495651391350007131;0.09608445831444172025;0.091400705452594849243;0.9808633719414890928;-2.3640155618252562952;-3.9683753071066560913;2.4716459238704158174;0.9103388722530324495;-1.6445206688436171394;1.7270617641242649309;1.3710565804733929607;-0.83839148097354254663;2.8957307008678130344];
            IW1_1 = [-2.4154308174831933265 1.1316734083225514773 -0.18798508993018392399 0.50752499754821289724 0.26643122967183774374 0.18496867771410097081 -0.084381354525935925448 0.07085802248326253383 -0.16205693380690422423 -0.77140106419434217866;-2.5214783316455919859 -0.61767172907865119935 -0.010085000823519651991 -0.17096230134548739965 -0.045858540570502259737 -0.075738671538680010786 0.011697215615493047544 0.89427052368563375584 0.017841955990787663339 0.36230068772579970826;-0.25273389398206069778 0.1998451513359776055 -0.012532406301889509326 -0.38710321317608464842 0.31744353291743737655 0.30392705893160959496 -0.48629036603932135341 -0.38527217222465176549 0.45637964870801228656 -0.0039633153238058721132;-0.14137397926041037066 0.22839219307467162334 0.69305240971952908335 -0.84866491694379209143 1.1527220931617192523 0.94439458980977164515 -0.30465769707891104945 0.37955421319014276405 0.40549188230617233542 -0.45339085921888461206;0.11430828401814195627 -0.11175304651558616575 -0.47704525296752076091 0.64561656700497560557 -0.916648058594788262 -0.7749657348705597526 0.31298468285909419873 0.54810711313383120302 -0.44114123533557358936 0.31640372701778679554;0.068730081107593654632 0.15546512469078524465 0.050894605481055565921 -0.18914112137188648921 0.018993525551347312241 0.19626799805238009933 -0.03532588208090964299 0.055668577385928771917 0.084229796875306933712 -0.26300932511266794656;2.6071271914854570184 1.4513736888071724351 0.0074584607403650983112 0.033216876158219783843 0.42244033655690699236 0.38798092469935352433 -0.15221425801286470048 0.38258859113757298642 0.73193461514120727873 -1.3507173453150382869;-0.028825152288339869061 -2.2081808752373786042 -0.061859892727623107256 -0.24630047120646228476 -0.14981246649827839601 -0.034852738621890928805 -0.030574791510504781278 0.037204871248042073462 -0.16536073768499864878 1.4525224141871446726;1.9848749283793558629 -0.65087570179592324493 -0.015401431195047120964 0.033833379986935939454 0.24766683536891354045 0.23619359560874905735 -0.061828332005150851702 1.3513681000858552839 0.474773459888372662 -0.10609490256059561641;-0.036487529426405596045 -0.97950882073050538068 0.13802254156185905787 0.51835177883852656677 0.5464992554107237499 -0.40436909977663293425 0.054607279876782710559 0.2782474473205402421 -0.15734213693163773273 0.19478347868992934577;1.0475016718182594833 2.6268066719483176286 0.10555368421258967682 0.39987028553426201549 0.18887696676464152401 0.070115438039725697106 0.059914126014487971428 -0.4542320291164952395 0.32303468224857806446 -1.6871729733077229785;2.7529274108906229834 -0.76501478843559456156 0.0018168068356873841133 -0.18797317063029567175 -0.056699649112672743934 -0.029790717506350518351 0.0098212456763821752437 -1.0800017164180586438 0.12325875968482226386 0.51340772757722163977;-4.5690876368247836936 0.29730258463782233136 0.021112838346440350457 -0.10558917971050130191 -0.12168972717260564953 0.0058857480734283767684 -0.081390986053868324968 -2.0464715164136473291 0.10785277117965144655 0.13108614636486201621;-4.8994890602218816866 0.35234830166711450516 -0.010013370059842452778 0.070144745130894969876 0.046229090500515740425 0.0041134165261851622294 -0.024505471792468629805 2.3295584707548151471 -0.18538945206563731127 -0.24224570259646180381;1.2615286870191546598 -0.11816400022351983279 -0.078271926206767192258 0.00040698798715488505471 0.12261647857164348352 -0.025018061970269584587 0.06275488348053152865 1.1711207094070865686 -0.15128987683747979753 0.011625717221515827204;0.25536351914056132362 1.5069046552312748144 0.040419707542304748882 -0.017552326786328720926 0.2518564757364556983 0.15937380698629941 -0.023214434568652780183 0.67048552164710273349 0.19329951559387478777 -1.2091523245907629391;1.4177390842659385317 -0.05695648275691699014 -0.10986006754139095165 -0.019963781991928088166 -0.13614685136084950234 -0.35015709512873327558 0.13999639673564417963 0.4251125746769744973 -0.8814843316401170803 0.1666711241041228253;0.41968287139809840047 -2.6683330419759339058 -0.098395196139792875933 -0.41495981268258286256 -0.11880388438818172137 0.011362744772546578068 -0.03973014734543630494 -0.50920599962215373768 -0.053892242227375410091 1.6413468033934477397;-1.3231026141444686139 -1.8493906165810243269 -0.056111174375088730681 -0.14865803431058338679 -0.14887176618909070402 -0.071243715652375777525 -0.020647144186653382941 0.35039232222030086694 -0.26197868120861045327 1.2119089090990931012;-4.9902268121632555875 -0.4172149111317444703 -0.0015895483878233280118 -0.15307047671214385476 -0.02357680456433861732 -0.066690710366615932325 -0.010463477723396974808 2.27970201198322453 -0.10341726132984332964 0.25867357616550795685];

            % Layer 2
            b2 = -1.3850185588017143168;
            LW2_1 = [-0.17239018609501124968 -2.4686570038789317216 -0.28132495741565127778 -0.65952833029270596654 -1.1372584731091543997 -1.9460586966551658428 -0.41574234276122307152 1.7184263126322130155 0.77527660655890684449 -0.38169495762906624492 -1.1742597402214667301 -2.7883027491604210901 -1.9903053748364492037 -2.069048211694330508 -2.037627117462171622 0.86637931156103764607 0.53376422641646292799 -0.94540883414461107659 -1.9486075433738543339 2.0027264449775534771];

            % Output 1
            y1_step1.ymin = -1;
            y1_step1.gain = 0.677942953516802;
            y1_step1.xoffset = -1.24456183373413;
            
            % ===== COMPUTING NN OUTPUTS ========
            xp1 = (regrN - x1_step1.xoffset).*x1_step1.gain + x1_step1.ymin;
            
            % Layer 1
            temp1 = b1 + IW1_1*xp1;
            a1 = 2 ./ (1 + exp(-2*temp1)) - 1;
            
            % Layer 2
            a2 = b2 + LW2_1*a1;
            
            % Output 1
            df1N = (a2 - y1_step1.ymin)/y1_step1.gain + y1_step1.xoffset;

        end
        
        % de-normalizing response (and also changing time units)
        df1 = (df1N'*norm.predScale + norm.predCenter)/(par.T);
        
    end
        
    % Form the DAE system
    diff = vertcat(df1);
    alg = vertcat(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26);

    % give parameter values
    alg = substitute(alg,p_res,par.p_res);
    alg = substitute(alg,p_s,par.p_s);
    alg = substitute(alg,T_a,par.T_a);
    alg = substitute(alg,T_w,par.T_w);
    alg = substitute(alg,T_r,par.T_r);

    diff = substitute(diff,p_res,par.p_res);
    diff = substitute(diff,T_w,par.T_w);

    % concatenate the differential and algebraic states
    x_var = vertcat(ER);
    z_var = vertcat(p_ai,p_wh,p_wi,p_bh,rho_ai,rho_m,w_iv,w_pc,w_pg,w_po,w_ro,w_rg,p_rh,rho_r,p_m,...
        w_pr,w_to,w_tg,m_ga,m_gt,m_ot,m_gr,m_or,V_p, mu_f,g1);
    p_var = vertcat(w_gl,mdot_p,GOR,PI,T);
    
end

