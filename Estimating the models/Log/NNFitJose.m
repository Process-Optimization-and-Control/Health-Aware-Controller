function [y1] = NNFitJose(x1)
%NNFIT neural network simulation function.
%
% Auto-generated by MATLAB, 12-Mar-2021 11:07:57.
% 
% [y1] = NNFit(x1) takes these arguments:
%   x = 10xQ matrix, input #1
% and returns:
%   y = 1xQ matrix, output #1
% where Q is the number of samples.

%#ok<*RPMT0>

% x1 = 1: sand; 2: pai_1_T, 3: pwh_1_T, 4: wro_1_T, 5: wrg_1_T, 6: prh_T, 7: pm_T, 8: wto_T, 9: wtg_T, 10: gLift_T;
% hiddenLayerSize = 20;

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-1.03797286522284;-4.46248802722717;-5.17852887458485;-4.48073715332605;-3.39270547151972;-4.89034267011056;-4.09638956423933;-22.1011864634595;-12.6725101583375;-1.68741224378933];
x1_step1.gain = [0.701732001985454;0.305647520874011;0.217139586607698;0.278754220747781;0.313094454989053;0.224952360464966;0.2398799497799;0.0883066374206199;0.130806471865733;0.582950679298246];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.6722792483485269255;-0.42794551761979493509;0.90823914038626818446;0.49728265501339979915;-0.54220161618642870316;-1.0574125697344736174;-0.18504744048835380932;-0.050635144090169392772;0.56812258991178576917;-2.2153704767782835816;0.26174524104815627412;0.30112132580966993745;1.1014891175005288027;0.028969840581056845175;-1.8729418612985893056;-0.69435213925219019249;-0.3837308175363531304;-1.1748024598933797247;0.24627934231708376034;0.63724273889490201661];
IW1_1 = [1.4425417287012727297 -3.9883829084529733855 0.010730207889524284379 -0.26087511669462959407 0.00071555176978392392127 0.1929206344952974217 -0.08457835754705340936 1.3882674613010050901 -0.11266323349720046287 2.1256799105111068293;0.53950726413098526457 0.98042095209695001579 -0.072326180367017448014 -0.19264948977040288747 -0.7299433338599803367 0.16234589514827740819 0.20536769637973770197 -0.16700914134727898319 -0.15387214119635495013 0.21170488067847903491;1.4883797442557591406 0.21379699735484336554 -0.67502021133950107501 0.42596822528159800081 -0.080413980624353134652 0.12089714842019801122 0.2444213940928378348 0.66331834771800313177 -0.35271582970743481367 0.047061116216261827949;0.27660994977060188171 -0.36333148955150279535 -0.015434539476298194186 -0.11773042664518294786 -0.23198610706058092612 -0.10738875474698106327 0.0008250566137152155928 -0.83528544755114331188 0.23908162344228486851 0.094878484703779561449;1.0939677295276322155 1.9582260166339275909 0.13483550720531492506 0.30122196302816262037 0.53127277306419051772 -0.050408570315821986385 -0.15593194521282896514 0.0845759271422918818 0.19120360593635615798 -1.6110528676613959931;0.82049261717892607582 0.64710984119649894719 0.023936832221339618409 0.11149750993057115644 0.11359672247777015397 -0.033961350931852443691 -0.017262525190507509859 -0.74200568656268106604 0.079774803185538409322 -0.51127280512098782772;-0.14588392083436976798 -0.82650725029636840091 -0.01616398430234436065 -0.1263574875752709048 -0.27875084277281780842 0.0081146685167722901189 0.025600227437710329492 0.67117352407405317383 0.1823668981504255604 0.27808877651170788026;-1.6543000019075964069 0.24211397437952347778 -0.78337839669612308047 1.0793415349672199621 0.15583152258488439101 -0.38776786060365453501 0.060168556897516538295 0.04880313886456840583 -1.3098622969093802748 -0.32263513477792715722;1.5787806193833684443 0.075863970020350779411 0.2994095845961399216 -0.54313286763907908661 -0.17762647743518455989 0.14845407991490994215 -0.011858167112330978907 -0.27980302544034757473 0.55706641135069212556 0.17962604094472972571;3.6293576309115231915 0.56453200961446192707 0.033745520568195318256 0.18663031482406389205 0.025884206613222802279 0.00076083097773667145679 0.052664684357108498824 -1.4053270458156399414 0.035806402723904812158 -0.32434737822438874932;-2.6819470595020180603 0.35908615924174780432 -0.10009327660579603325 -0.085503662620257853133 0.21002203387358417408 0.067321600699544764668 -0.10377396641765268548 0.41716832491653443871 -0.050910293393761650582 -0.34735299708086464499;-2.1400108264990520546 -0.36890850083836712603 -0.04061067971394231868 -0.14420972347364929833 0.029278839559895974282 0.020689502357154793827 -0.063122865671661976261 0.60664389654316586142 0.024195708112017207808 0.15894673116526153445;1.1813960117875423883 0.30796413994086557775 -0.81956146191702627934 0.65653870813802583939 -0.38092038556262614035 -0.38523896302273646919 0.094687142714984767755 0.21995320828529155155 -1.0703663241817567009 0.14648301769787225313;0.22460450899689690663 -0.87259378125932107739 -0.13377108587569772058 -0.20664321614221772827 -0.53738256862694933069 -0.019944294478215572608 0.099774879818675443754 -0.091046492415845448054 0.36329849380191708841 0.32547688301314725567;3.252669028343739388 -0.61597473698053817337 0.053095862690751338686 -0.0074618773303536207672 -0.10985089644490575456 -0.050317206288186901164 0.05108995998023423829 -0.83267601077554143085 0.14533441026414894992 0.39999955105175483139;-2.3581990198703608463 -0.86076923513364567153 -0.050473574688998691651 0.32824956624567053209 0.069667431685819750098 -0.025812355306680301315 0.017051961717360807069 0.21295450881853453251 -0.17447205819835051588 0.28866738635094102428;1.1386693769988711633 1.0182385601521166674 0.1790947455579601777 0.41018547870371746544 1.033705601329207191 -0.23365590964042448618 -0.3083207323853972559 -0.037051847330360751509 0.36336636923046322645 -1.6098230721156510103;-0.70924097126931728852 0.72202656743417537033 0.06995341236889505987 0.10029341507785689835 0.30679086649469133885 -0.055036431553584348741 -0.068778987122707135882 -0.3983535936478306283 -0.16772380436014538985 -0.2412839801337656187;-0.98200668092497944706 -2.7720463429816617484 -0.068318478758515405191 -0.27091837302357946182 -0.18069723469737877797 0.033353963235018130318 -0.0069788104583871843345 0.18089543594817231176 -0.18078210716272824499 1.7190610573794788962;-0.89943544630851235411 -0.99895656182084424657 0.039905542902102374148 0.019986489047563813642 0.60798262716323692167 -0.12662762691075316046 -0.14713008573380620403 0.48690525237197040953 0.14917100939088104306 -0.02035097048936793801];

% Layer 2
b2 = 0.15800871537246910514;
LW2_1 = [0.23886473571332755839 -0.62896788593465646766 -0.26600674569424709315 -1.0719876903795724221 -0.77611083951924308977 1.8612575279113681148 -1.2882285871612131345 0.51038022634394430987 1.8604005924431143981 -1.6107154171239284057 -1.219981394823503873 2.2671897628579182893 0.27835206505983323133 0.52211092641970247019 1.5584271883574194462 0.62686941298879428341 0.35363753555170018839 -1.2975448112979310888 -0.64266758907385779231 -1.0857898056481238847];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.677942953516802;
y1_step1.xoffset = -1.24456183373413;

% ===== SIMULATION ========

% Dimensions
Q = size(x1,2); % samples

% Input 1
xp1 = mapminmax_apply(x1,x1_step1);
% computes: xp1 = (x1 - x1_step1.xoffset)*x1_step1.gain + x1_step1.ymin

% Layer 1
a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*xp1);

% Layer 2
a2 = repmat(b2,1,Q) + LW2_1*a1;

% Output 1
y1 = mapminmax_reverse(a2,y1_step1);
% computes: xp1 = (a2 - y1_step1.ymin)/y1_step1.gain + y1_step1.xoffset
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end