function erosionHat = Diagnostics(xPlant,erosionHat_1,yMeas,uk,sandMeas,par,modelFlag)
% Diagnostics layer to estimate the current value of erosion

    %% normalization (for stepwise and NN)
    norm = load('normalizationValues');
%    load('StepwiseLinear');

%%
 if modelFlag(1) == 1
        % Use phenomenological model - perfect information
        erosionHat = xPlant;
    else
        % computing the regressors
        %regr = [mdot_p; p_ai; p_wh; w_ro; w_rg; p_rh; p_m; w_to; w_tg; w_gl];
        % separating the regressor for the three wells. Some regressors are
        % shared by the three models
        regr = [yMeas(1), yMeas(2), yMeas(3); %pai 
                yMeas(4), yMeas(5), yMeas(6); %pwh
                yMeas(7), yMeas(8), yMeas(9); %wro 
                yMeas(10),yMeas(11),yMeas(12); %wrg 
                yMeas(13),yMeas(13),yMeas(13); %prh
                yMeas(14),yMeas(14),yMeas(14); %pm
                yMeas(15),yMeas(15),yMeas(15); %wto 
                yMeas(16),yMeas(16),yMeas(16); %wtg 
                uk(1),uk(2),uk(3)]; %wgl        

        % normalizing regressors with data from (ExpDataAnalysis.m)
        regrN = (regr - norm.regrCenter)./norm.regrScale;
        
        if modelFlag(2) == 1
            % number of regressors
            nReg = length(regr);
            
            % modeling stepwise 
            % putting the weights in the right place
            SW_M1 = zeros(1,nReg);
            SW_M1(1) =  0.0755381434097584;
            SW_M1(2) =  0.000243142718240501;
            SW_M1(3) =  0.0100252738194035;
            SW_M1(4) =  0.00380047454552493;
            SW_M1(5) = -0.000143628353752839;
            SW_M1(7) = -0.0114289976747677;
            SW_M1(8) = -0.0298136949369108;
            SW_M1(9) = 0.927543516164127;
            
            SW_M2 = zeros(nReg,nReg);
            SW_M2(1,3) = -0.0125670258873362;
            SW_M2(1,5) =  0.0151494785810274;
            SW_M2(1,7) =  0.0128678467101190;
            SW_M2(1,8) = -0.0153404875401273;
            SW_M2(1,9) =  0.0138351276033684;
            SW_M2(2,5) =  0.00109630362172230;
            SW_M2(3,9) =  0.0117344652682547;
            SW_M2(5,9) = -0.0154268998662768;
            SW_M2(7,8) = -0.00266142658066635;
            SW_M2(7,9) = -0.0101639466718962;
            SW_M2(8,9) =  0.00898894223338620;
            
            intercept = -0.00919521626422190;   

            %building model
            erosionRateHatN = [];
%            erosionRateHatTest = [];
            for well = 1:3
%                erosionRateHatTest = [erosionRateHatTest, StepwiseModel.predictFcn(regrN(:,well)')];  
                erosionRateHatN = [erosionRateHatN, intercept + SW_M1*regrN(:,well) + regrN(:,well)'*SW_M2*regrN(:,well)];
            end
 
        elseif modelFlag(3) == 1
            % Use neural net model
            % Input 1
            x1_step1.xoffset = [-4.46248802722717;-5.17852887458485;-4.48073715332605;-3.39270547151972;-4.89034267011057;-4.09638956423933;-22.1011864634595;-12.6725101583375;-1.68741224378933];
            x1_step1.gain = [0.305647520874011;0.217139586607698;0.278754220747781;0.313094454989055;0.224952360464966;0.2398799497799;0.0883066374206199;0.130806471865733;0.582950679298246];
            x1_step1.ymin = -1;
            
            % Layer 1
            b1 = [2.1798274445631249385;-2.0573794874510471509;1.8615123586554185309;-2.049260576459216221;-1.8629720127482731762;1.8000929933905742697;-1.8263040091667326337;-1.4623553313527706266;1.4565479669777872118;1.2699644299198349362;1.5815778994964344584;1.2025913624444877215;0.98560019584480118571;-1.0628422847694944942;-0.98161329145620701553;0.82722121399777925888;0.78736841272189916907;-0.66481747429254123283;0.38418470249910974257;-0.33386602403848009146;0.62286482332406001294;0.43411831369831255234;0.13634563146247311027;-0.00055379147205731043407;-0.060237286145397542225;-0.21928329327698811713;-0.16464372547631736521;0.12005641044072273838;0.28419191153564138386;-0.56213023013270135841;-0.46396233979125084401;0.91440446524627594371;0.88236962343122893326;-0.78843004563356411385;-0.83338286085928581226;0.85162889512153749916;1.0983098003887197258;-0.84207265278513865336;-1.0752921997979758562;1.1764902667746457432;-1.3747461615091365328;-1.7519107369151605269;1.5674989748981424942;-1.6080916773306022094;1.9021536144518114142;1.955728199228989439;2.3207579197318932707;2.2968329250327803592;1.9061899849765546744;-2.1477467307203399471];
            IW1_1 = [-1.1590386666724961096 -0.79699349001173525053 -0.77248864963340702072 -0.71217593624325103985 0.67641786306272211338 -0.54800415716663919419 -0.18409046951114155233 -0.14832288918795283328 -0.63113787591284109535;0.098797229937017325141 -0.81488305312285758131 0.62279338852340448973 0.98130851337816760172 0.70140605512246867459 -0.6573650208793230254 1.0750077709405356785 0.30233244722029573026 0.0096182402659141658097;-0.17635659706743611608 -0.58136843394441251043 1.1232222854868449957 0.53101706462729036939 0.85246908448098257516 0.32272031662518380069 -0.73645314690601015872 -0.52115087872079524089 0.28856275789931140618;0.49225832992995421922 -0.3994821335845434862 0.65648045093042073361 -0.27899861332624331656 0.005039269485490011774 0.022317025709633085684 0.85135341095263605826 -1.4086189070559858028 -1.0456466843118505938;0.19016268274677639027 0.84022534532895709614 1.0862395905043114031 0.81307428162702743979 -0.7390605004795952615 -0.903246531213847903 0.29722299805743168655 -0.23909660874931196184 -0.50776411190367753967;-1.0292020699193995537 -1.0284240857677988679 -0.72323803027460864534 0.37208853414762299572 -0.59475791652411724719 0.4221580816587881535 -0.88514557491291079305 -0.37241012271870244099 -0.099580672416987758844;0.20246151249354829571 -0.32296738159542914826 1.1797692588289763194 0.60648960096269288744 -0.79029833925059844724 -1.1042945328256281723 0.23232873695813446857 0.059281285113776636952 -0.96299333097251338387;0.51268570992100981343 0.45183825666325350134 0.25647405227440744913 0.31026653657353719939 -0.1362179235796261989 0.18218649205857051498 0.97050403673971497032 1.0104645141247559081 0.62559324736000321288;-0.81159293912694485673 1.0428801008072525747 -0.31299708012949539748 -0.7452545787618123363 0.44917031001882012919 -0.94520331125958190199 -0.27400158142352765145 -0.43863069480757677088 -0.65688533400532378437;-0.56258008349708943374 -1.0017970249375542924 0.18527065261407191232 0.65212918825388033195 -0.70440040505697576467 0.43083973535443609304 0.39182674089317265365 0.11496002386879981827 -1.004010625580097793;-0.090683472038655787983 -0.60545054452687085966 0.5803037981735260642 0.52839542210652346022 0.60119081681072483736 0.67469186064477126852 1.0011872659222740278 0.94867152509506580405 0.54484016891390940618;-0.73336317550733931547 0.70657009118205182929 -0.72573944496047471109 0.23003775785568816858 -1.1279122661328948585 -0.0079771365055598380395 -0.3619180587553524342 1.2210240752885037629 0.27080105913943880447;-0.13508071842691662945 0.58749143030319972425 0.794270855194250891 0.97603607764239996847 1.015779433542914223 -0.37448745101425084947 0.67988053976004603207 -0.27307705937745019309 0.87636221984866824464;0.30746355056318419585 -0.035697213509789285646 0.31360561570836831269 0.6242206398246903376 0.033334443081142556276 0.43768245071047023886 -0.62871041731231447347 0.89706555783432528894 0.21775858193171759858;0.57818243866953011878 0.092952369660553743613 1.0808308374326520696 -0.41815470514646724265 0.37300783817356580885 0.52956348727930191433 0.98133067677475316337 1.0822327390784785006 -0.58263385061369976636;-0.15759849128533379159 -1.0116060899937944217 -1.1463858012686147791 -0.99133745934940353806 -0.72745118484255966163 0.52277399092078913334 0.53015949481311608338 -0.024324509749098666145 -0.23623937755144469586;-0.1540885236917124701 -0.5978948040923621754 0.23603355399941655657 0.58129951248796207786 0.83867595079379886336 -0.31114578162429645802 -1.2106349709387185243 -0.067207336139262704999 -1.3598663036447100616;0.47409127780329168855 -1.1728700974120278566 -1.187022471096889964 -0.93762365797092306341 0.060647467978965030611 0.23905694210418718115 0.56294996041728062597 -0.46089421618771014844 -0.087074892047195845413;-0.46725444095241214226 -0.74082475780320378167 0.40402922430110116414 1.0251041106393172786 -0.78686162954815230819 0.29427209655605313099 0.056125883980321089606 1.0587485609597910852 -1.0732814622423079154;0.4456154596460165096 0.11869961166749141979 -0.85708101986257212435 -0.35530739718530379712 -0.84397015384473306199 -0.23267432661140749484 -0.028255872614553041722 0.88157598352422716914 -0.6942708272141653314;-0.26053131549916203058 -0.028319020203103305311 -0.26539672452994483853 -0.23348474620382930755 -0.46481159257776583082 -0.16532320207394390899 -0.13111830213911515597 -0.21857280335606038557 0.68640874088595338343;-0.54245489522821932482 -0.50276578037344066274 0.6760617078531055979 -0.81896795093647400954 0.60636990993012174211 -0.68990448772834911395 -0.57655779256111883146 0.56106983453155767627 0.21719811137118577582;-0.099592574511474116217 -1.3326564361991268814 -1.0008405691408259752 -0.058318691555542992444 0.42370414233587189168 0.44654224221809096784 -0.30202671055388030652 -1.1576996690960668079 -0.03080188324094252153;0.081163005686760394952 1.3758887457715249791 0.4795387998833334775 0.056949174646781244857 -0.21513358554156841729 0.16390561258040797776 -0.52230220310806596551 0.6936890175094312383 1.1930565664416052929;0.4501697684759001783 0.67316331516504424037 -1.0047334004729275847 -0.57122922802725883962 -0.48108883653252781931 0.37660454398994125391 -0.087450010780991196957 0.11882647063421380995 0.84619609823674690752;0.90570376207389835521 -1.02882117213961366 0.41413446144262394677 -0.37854805696946658422 0.51962854565924654349 -0.43302618469250542255 -1.0545694901872886273 -0.85990005785007872507 -0.28063378326706556187;-0.26864011613968657999 1.2232365486545464517 -0.16126931997324706058 0.76249662153164521339 0.30737451478160021656 -0.25976599434303609959 -1.024253778132076631 1.0398996355427949645 -0.46680849521784822942;0.72233378072553422466 0.65742315652676186488 -0.79537268179348930008 0.5997823550155212402 0.65753830673976587651 -0.67033865664737113388 -0.30897174657443515766 0.01746467871045798495 0.57566341846638668933;0.26272682168785643775 -0.90059789657023370246 -0.98233957412910399842 -0.47020110469776693618 1.3101781128324858638 0.23891763552452316421 -0.55421213982746431181 0.0012085939875530230751 -0.55184573393443092648;-0.48077489572139708862 -0.98351444032327262068 0.4848421948156674266 -0.41603663427100529004 0.9303958852055345119 -0.77582085163798131244 -0.98376382712210408421 0.2046746302419850394 0.27807650596349919869;-1.2011476824257718565 0.97005140229132180707 -0.26401196076628069553 0.13152647152278995613 0.53613313708781218381 1.3658348726506543969 -0.17107298891360062099 0.048455327266491696192 -0.058916655140559402726;0.59964245157308315726 -0.10367648163528757344 0.64206457835769703557 -0.68020311930778643017 0.81550160107737523774 0.63285251206116321931 0.67917887331007087415 1.0621127085442878535 -0.3242293328644266226;0.26068952267809941681 -0.32296246434087882982 0.10029302793129915572 -0.77252203360896509743 -0.019996875778740326457 0.63053844247206125573 1.1843801915047771356 1.3152776654653934862 0.504252208325274065;-1.0070400909576369131 0.29293445567421505515 0.59790745254008959009 0.90256952957622271771 0.73445213635008443909 -0.098216045095856407965 -0.82986508795633051605 -0.46824755667685269422 0.39215556269885670782;-0.95379691101366048667 -0.6883891940167005119 -1.0685577246427211406 -0.040298245933887422154 0.65763722228855714391 0.17276876007922462031 0.49462053171365705229 -0.84277927620841985856 -0.82651162556128554559;0.30814557786734853018 -0.58261071447689605218 -1.3915534030866831827 -0.12432239821110645239 -0.88035497608871471353 -0.596879548491172085 0.65669956185986055885 0.51802486864884877971 -0.54826775095672941607;0.8631998640808919987 -0.93450964693359339641 -0.18072238081074881721 -0.53820044651640808198 -0.41631960908563225843 0.93276967504066932158 -0.64686380820280897641 0.71486602573181456943 0.68709920771537080775;-0.58936579628755081384 -1.1768071816518543837 -1.244018538866214163 -0.502791159761638351 1.1517262278749618165 -0.45025220285367650863 -0.097916181235936108895 0.0039371071567724479079 0.49359195577059261639;-0.80431524585474700739 1.4322166283071871007 0.75523692298971012438 0.59054477775606384338 -0.64077346303498994295 0.0088870348969691120855 0.56998535269912842249 -0.56311479027823219123 0.19000608324756124157;1.0381150021960867669 0.9546819436980813478 1.0183971109478029415 -0.24196286126545790696 1.0280212979329439449 -0.38296441911574080619 -0.14989156087946542417 0.28634871755540142191 0.42176036772014496856;-0.54545550867830350228 -0.69971776083748438158 1.1649367548491393887 0.062791969524111687639 0.73463144755474241165 -0.81726081890696922194 -0.82235925508366647385 0.56443080497383957272 -0.093207494908287516489;-0.016126077805874713261 0.10944553675823007199 -0.76045641327351232164 0.47315008008830117259 -0.47031522288295884993 -0.56857610873547592689 0.087381821872159079478 -0.49314643328061807948 1.0175303162633257958;0.13360738970949584403 -1.0103731233704480452 -0.18588648253842002256 0.20184291971541298838 0.89784452609993792294 -0.43081673241041223976 -0.98223587740829643167 -0.5638851710863774791 -0.067319491827354152602;-0.95167564970700102567 0.25081598941185506702 0.82829256594133138236 0.28321169318609834553 -1.009733995817717167 -1.2489738023096454977 0.70810618293345806151 -0.60389040435191210232 0.19401433343130608522;-0.54149005750730128117 -0.94660398896892183895 -0.5364021493647652683 -0.90936871695515220093 0.56540838733706144659 0.70795515097368855173 -0.38939499378821063047 -0.0033874140447711295621 1.4492219537785078209;0.84204887322754917722 0.86849312263024491809 -0.76523587041472185266 -0.066552958189353084162 0.80204235678818180943 -0.41481389173226224099 -0.2376954602465767008 0.44134771840023778733 -0.95178964262203835656;0.76367475648548854839 0.84746708163497852517 -0.43161605835272415232 -0.69043946027208535909 0.79687466263100259489 -0.6548387386362722884 -0.025468863761413173402 0.79335589030070641225 -0.53201302204178602118;0.72173704558330620262 0.83451066159322218052 0.55550616212741343869 0.42954169282158211063 -0.57850840433869787915 -0.95618396782104930764 -0.52892270389017725929 1.0765454038965949213 -1.0970026244098884938;1.2356392959677684473 1.1688526478948126908 -0.58122789031638288915 -0.46834435512901195775 0.63081560836912409229 0.64734891337441069492 -0.29476891708776897794 0.79282386526535131299 -0.19809493214193973909;-1.0135933448356300168 0.14339307680963195302 0.96711990983754436702 0.027037496254616384672 -1.0538476898804203241 -0.24235106847597082114 -0.52909361921150732666 -0.97303197284450615445 -0.23484293128128916339];
            
            % Layer 2
            b2 = 0.30012127919607584081;
            LW2_1 = [0.0128790468769139025 0.025675922910405296617 0.24233096859867975237 -0.33241959376440294305 -0.12506139117453857201 -0.14756586316505460044 -0.2514343499159332751 0.13911101490080590715 -0.013798477187762800064 -0.29816878361017207233 0.16999313202454069405 -0.13055635527823777897 0.035770761579104796979 0.21471641997370161103 0.13346107209388907711 -0.014699422914239873225 0.025414357650865519733 -0.0087784485239775656568 0.16996138437626368001 -0.13611627299089265608 0.96523676675784564338 -0.070042475938084824971 -0.085903171421176247047 -0.023829872320772729766 0.095462721866351907263 0.15224447622567510718 -0.094738817110316023151 0.06172767363915163813 0.055084245084018303162 -0.4995770666799899673 0.027323222861291029256 0.19462378897993298787 0.20516088646458338629 0.36362443426605067787 0.005193521025020380287 -0.078738520307218631822 0.086945197893806017153 0.15739518353443754384 -0.055230333144792864708 0.081002757000515360208 0.49359597924561149362 0.33305627983936381131 0.071598162921745436837 0.17683952168754082934 -0.2242694865749460531 -0.37748839411890205975 0.61802347987722972 -0.58483574461309129067 -0.86062218590467531865 -0.70965992896541552071];
            
            % Output 1
            y1_step1.ymin = -1;
            y1_step1.gain = 0.574079476354099;
            y1_step1.xoffset = -1.68347654001744;
            
            % ===== COMPUTING NN OUTPUTS ========
            xp1 = (regrN - x1_step1.xoffset).*x1_step1.gain + x1_step1.ymin;
            
            % Layer 1
            temp1 = b1 + IW1_1*xp1;
            a1 = 2 ./ (1 + exp(-2*temp1)) - 1;
            
            % Layer 2
            a2 = b2 + LW2_1*a1;
            
            % Output 1
            erosionRateHatN = (a2 - y1_step1.ymin)/y1_step1.gain + y1_step1.xoffset;

        end
        
        % de-normalizing response (and also changing time units)
        erosionRateHat = (erosionRateHatN'*norm.predScale + norm.predCenter);%;/(par.T) --> don't need to adapt the units
        
        % erosion level estimation
        erosionHat = erosionHat_1 + erosionRateHat;

        % snipping values smaller than zero  
        erosionHat( erosionHat < 0 ) = 0;

    end

end

