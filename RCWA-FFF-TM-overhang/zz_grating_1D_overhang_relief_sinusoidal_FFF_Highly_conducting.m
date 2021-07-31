% planar diffraction by a sinusoidal relif grating
clear all
%clc
lam    = 1;
%permW = 1.0;
%permW = 1.72^2;
%permW = 3^2;
%permW  = (1+5i)^2; 
permW  = (0.1+10i)^2;   %NOTE that the permitivity is complex!!!
%permW  = (1+10i)^2;   %NOTE that the permitivity is complex!!!
thI    = 15*(pi/180);%+1e-3*(pi/180);
epsB   = 1;
Lam    = 425e-3/525e-3;
d      = 360e-3/525e-3;
epsW   = permW;
epsS   = epsW;
k0     = 2*pi/lam;
q      = lam/Lam;
q0     = sin(thI);

%Phi=0*(pi/180);
Phi=30*(pi/180); 

%% Control parameters
StrucParam.CS='CC_CS'; 
                                      %'CC_CS': default,cos^2 and cos*sin: according to Popov, et.al. , and Habib Mohamad et. al.
                                      %'C_S': cos and sin: according to T.V.
                                                                          
StrucParam.filtering=1; %0: default, no filtering; 1: filtering; According to Nikolay M. Lyndin, et. al.                     
StrucParam.threshold=6.5;  %for eigenvalue
%Near field plot;
StrucParam.PlotField='NO';  %'NO':default,not plot field;
                                        %'YES',plot near field, quite slow.
StrucParam.h_upper=1;      %h_upper=StrucParam.h_upper*h_grating;
StrucParam.h_lower=1;      %h_lower=StrucParam.h_lower*h_grating;
StrucParam.Number_of_Period=2;
StrucParam.Resolution_x=200;
StrucParam.Resolution_z=200;                                        
%% Functions of sin and cos tangent to the profile

K1=2*pi/(Lam*cos(Phi));

%%{
%NVM
%simple profile which can be discribed by a single-value function
%exaple: simple sinusodial profile.
%fsin=@(x) 1./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce cos(phi(x))
%fcos=@(x) (((pi*d)/(Lam)).*(sin(K1.*x)))./(sqrt(1+((pi^2*d^2)/(Lam^2)).*(sin(K1.*x).^2))); %predpis funkce sin(phi(x))

%comple profile described by parametrized function, like circle, overhang sinusodial function etc.
%Parametrized derivatives.
%fsin=(dx/dt)/norm;
%fcos=(dy/dt)/norm=sqrt(1-fsin^2);

if strcmp(StrucParam.CS,'C_S')

fsin=@(t) (sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t)))./sqrt((sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2+(cos(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2); %predpis funkce cos(phi(x))
fcos=@(t) (cos(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t)))./sqrt((sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2+(cos(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2); %predpis funkce sin(phi(x))

elseif strcmp(StrucParam.CS,'CC_CS')
    
%NOTE that here cos and sin have been exchanged relative to the definition above
%fcos=C^2; fsin=C*S
fcos=@(t) (sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2./((sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2+(cos(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2); %predpis funkce cos(phi(x))
fsin=@(t) (sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).*(cos(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t)))./((sec(Phi)+sin(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2+(cos(Phi)*(K1*(d/cos(Phi))/2).*(sin(K1.*t))).^2); %predpis funkce cos(phi(x))
else
    disp('The definition of cos(x) and sin(x) is wrong!');
end
%}

%{
if strcmp(StrucParam.CS,'C_S')
% No NVM factorization, i.e. just the norm staircase approximation with Li's rule.
fsin=@(x) 0.*x;
fcos=@(x) 1+0.*x;
elseif strcmp(StrucParam.CS,'CC_CS')    
%NOTE that here cos and sin have been exchanged relative to the definition above
%fcos=C^2; fsin=C*S    
fcos=@(x) 0.*x;
fsin=@(x) 0.*x;
else
    disp('The definition of cos(x) and sin(x) is wrong!');
end
%}

%% Preallocate fields for reflection amplitudes and errors
errorS_R=zeros(1,1);
errorP_R=zeros(1,1);
errorS_T=zeros(1,1);
errorP_T=zeros(1,1);

%% Reference values for sinusiodal grating, computed by C-Method
%%n=1+5i
%RP_ref=[0.0674606805767746 0.00319796649824055 0.303105896615171 0.124877473586179];
%RS_ref=[0.327219701682201 0.207465501464766 0.117286155151613 0.0792741802547489];
%TP_ref=zeros(1,6);
%TS_ref=zeros(1,6);

%{
%%n=1.72,theta=0 [deg], Phi=30 [deg];
RP_ref=[0.001297134015254];
RS_ref=[0.010174769316324];
TP_ref=[0.312353806742803   0.681667737517812   0.004681364089045];
TS_ref=[0.726658666521634	0.229551910187410	0.033614676815176];
%}

%{
%%n=1.72,theta=15[deg],Phi=30 [deg]; C method,WRONG,sort beta, need rotation
RP_ref=[0.004048697696454   0.002273649038408];
RS_ref=[0.020351069199529   0.002595456005422];
TP_ref=[0.204563103140346   0.715328171657040   0.073786421337080];
TS_ref=[0.465951079541316   0.429741994958952   0.081360436389336];
%}

%{
%%n=1.72,theta=15[deg],Phi=30 [deg]; LPEM
RP_ref=[0.013847619286624   0.001321278113111];
RS_ref=[0.033479360118408   0.003225124752801];
TP_ref=[0.276865218496503   0.698838843903608   0.009143687581078];
TS_ref=[0.429215218389913   0.503996524935058   0.030083766291384];
%}

%%{
%%It seems that this grating is not lossy enough, so maybe we need to?????
%%give n2 a larger imaginary part test the different algorythms.
%%n=1+5i,theta=15[deg],Phi=30 [deg]; 
%FDTD %PLEASE DONOT DO THIS TEST.
%RP_ref=[];
%RS_ref=[];
%LPEM
RP_ref=[0.690696929434143   0.044679728260437];
RS_ref=[0.271923709974736   0.502654339036197];
%FEM
RP_ref=[0.6969397576849351	0.042593465794829075];
RS_ref=[0.27262164737017386	0.5021319693110445];
%Stardard FFM:RCWA-1D+S method:N=401,M=200;
RP_ref=[0.685220257896775 0.0451030975250511];
RS_ref=[0.271998824578871 0.502563618626655];
%Stardard FFM:Reticolo+modified T method:N=401,M=200;
RP_ref=[0.685220257895786 0.045103097525231];
RS_ref=[0.271998824578443 0.502563618627200];
%FFM-FFF:N=401,M=200;
RP_ref=[0.685220257896644 0.045103097524921];  %Staircase_CS

RP_ref=[0.685220257896636  0.045103097524920];  %Staircase_CCCS

RP_ref=[0.679735717417338 0.043012742570048];  %FFF_CS,layer by layer

RP_ref=[0.680169300567980 0.043733149380730];  %FFF4_CCCS,layer by layer


%FFM-FFF-new:N=61,M=100;m=12;
PR_ref=[0.593521420734058 0.033795070317326];  %Standard RCWA-1D+S
RP_ref=[0.593521420734065 0.033795070317335];  %Staircase

RP_ref=[0.583566274577521 0.032372514988748];  %FFF,once

RP_ref=[0.667510806093971, 0.034303475841890];  %FFF_CS,layer by layer

RP_ref=[0.668736568932016 0.041550849121362];  %FFF4_CCCS,layer by layer

%FFM-FFF-new:N=61,M=200;m=12;
RP_ref=[0.580165236582215 0.033063285468699];  %Staircase

RP_ref=[0.681570672315945 0.036551832260432];  %FFF_CS,layer by layer

RP_ref=[0.682235478195283 0.042372902878796];  %FFF4_CCCS,layer by layer

%FFM-FFF-new:N=61,M=250;m=12;
RP_ref=[0.576906391756812 0.032716776594237];  %Staircase

RP_ref=[0.684487204014025 0.036991332298476];  %FFF_CS,layer by layer

RP_ref=[0.684576397729647 0.043444529012805];  %FFF4_CCCS,layer by layer

%FFM-FFF-new:N=61,M=300;m=12;
RP_ref=[0.577917632401220 0.032938555643711];  %Staircase

RP_ref=[0.685968647391606 0.037153986395398];  %FFF_CS,layer by layer

RP_ref=[0.686557503547818 0.043808361996801];  %FFF4_CCCS,layer by layer

%FFM-FFF-new:N=81,M=200;m=12;
RP_ref=[0.611794654180626 0.036396764482898];  %Standard RCWA-1D+S
RP_ref=[0.611794654180615 0.036396764482898];  %Staircase_CS
RP_ref=[0.611794654180613 0.036396764482897];  %Staircase_CCCS

RP_ref=[0.608290185707760 0.035151844834948];  %FFF,once
RP_ref=[0.682107033399329, 0.038013455352558];  %FFF_CS,layer by layer
RP_ref=[0.681599347907159, 0.042963578909171];  %FFF4_CCCS,layer by layer

%FFM-FFF-new:N=101,M=300;m=12;
RP_ref=[0.626801354406051 0.037930934014082];  %Standard RCWA-1D+S,W-S
RP_ref=[0.626801354406098 0.037930934014083];  %Staircase,C-S,W-S


RP_ref=[0.624271398463843 0.037365482267794];  %FFF,once

RP_ref=[0.684959865382469 0.039917326428615];  %FFF_CS,layer by layer
RP_ref=[0.685056825574264 0.043589805451875];  %FFF4_CCCS,layer by layer

%%FFM-FFF-new:N=101,M=200;
RP_ref=[0.629234899653677 0.037974148151751];  %Standard RCWA-1D+S,W-S
RP_ref=[0.629234899653715 0.037974148151734];  %Staircase,C-S,W-S
RP_ref=[0.629234899653714 0.037974148151734];  %Staircase,CC-CSW-S

RP_ref=[0.629234899653677 0.037974148151751];  %Standard RCWA-1D+S,W-t-S
RP_ref=[0.629234899653721 0.037974148151735];  %Staircase,C-S,W-t-S
RP_ref=[0.629234899653724 0.037974148151734];  %Staircase,CC-CS,W-t-S

RP_ref=[0.680447861919047  0.039006330177836];  %FFF_CS,layer by layer,W-S
RP_ref=[0.678909180209482 0.042473568386308];  %FFF4_CCCS,layer by layer,W-S

RP_ref=[0.680447861919065 0.039006330177838];  %FFF_CS,layer by layer,W-t-S
RP_ref=[0.678909180209494 0.042473568386311];  %FFF4_CCCS,layer by layer,W-t-S

%C method
RP_ref=[0.691131155045779   0.044826802871665];
RS_ref=[0.271937696305567   0.502697416029626];
%}

RP_ref=[0.956070372866905 0.036126429767618];
%RP_ref=[0.897562186983547 0.028743965469486];

%{
%%n=1.5,Phi=0 [deg];
RP_ref=[0.002718305253944];
RS_ref=[0.020885635586605];
TP_ref=[0.158865089352025   0.679551452063414   0.158865078714946];
TS_ref=[ 0.080711339427445   0.817691514993539   0.080711337946249];
%}

%{
%FFM-FFF-new:N=61,M=100;m=12;
RP_ref=[0.247361277798774 0.199847804966880];  %Standard RCWA-1D+S
RP_ref=[0.247361277798784 0.199847804966873];  %Staircase

RP_ref=[0.312485873264685 0.293510393773617];  %FFF_once
RP_ref=[0.312485873264685 0.293510393773617];  %FFF_CC_CS
RP_ref=[0.305372515504232 0.292439782054019];  %FFF_CC_CS

%FFM-FFF-new:N=61,M=200;m=12;
RP_ref=[0.244308773768562 0.196364571403490];  %Staircase

RP_ref=[0.319556706473012 0.300399457139619];  %FFF

%FFM-FFF-new:N=61,M=250;m=12;
RP_ref=[0.242554087576985 0.197907809658782];  %Staircase

RP_ref=[0.321661078820742 0.302300006185708];  %FFF

%FFM-FFF-new:N=61,M=300;m=12;
RP_ref=[0.243652848768707 0.197698705219392];  %Staircase

RP_ref=[0.322326204987162 0.302964795714563];  %FFF

%FFM-FFF-new:N=81,M=200;m=12;
RP_ref=[0.261151456104021 0.220954404598872];  %Standard RCWA-1D+S
RP_ref=[0.261151456104022 0.220954404598859];  %Staircase

RP_ref=[0.318055645643127 0.300508461674967];  %FFF_CS,once
5
%FFM-FFF-new:N=101,M=300;m=12;
RP_ref=[0.268926202563466 0.231562690969994];  %Standard RCWA-1D+S
RP_ref=[0.268926207190528 0.231562695266635];  %Staircase

RP_ref=[0.319618184133119 0.302145013133523];  %FFF

%FFM-FFF-new:N=201,M=100;m=12;
RP_ref=[0.305912743287923 0.281613093473929];  %Standard RCWA-1D+S
RP_ref=[0.305912741321618 0.281613107429333];  %Staircase, precision convert
RP_ref=[0.305912743286996 0.281613093474450];  %Staircase, NO precision convert

RP_ref=[0.300006706643219 0.278596604689877];  %FFF

%FFM-FFF-new:N=201,M=200;m=12;
RP_ref=[0.311495606156169 0.292450232286330];  %staricase
RP_ref=[0.308966344404275 0.278448183500382];  %staricase

RP_ref=[0.311495606156169 0.292450232286330];  %FFF,once
RP_ref=[0.311495606156169 0.292450232286330];  %FFF,layer by layer

%FFM-FFF-new3:N=201,M=200;m=12;
RP_ref=[0.325722911382555 0.054356829805529];  %FFF

%FFM-FFF-old:N=401,M=200;m=12; By the author.
RP_ref=[0.313336500143800 0.290017147531323];  %FFF
RP_ref=[0.313055826576037 0.289669289721489];  %FFF2

%FFM-FFF-new:N=401,M=200;m=12;
RP_ref=[0.314471651780726 0.290762212267354];  %Standard RCWA-1D+S
RP_ref=[0.314471649364373 0.290762219552851];  %Staircase
RP_ref=[0.312745729594486 0.292558328347222];  %FFF

%%n=1+5i,theta=15[deg],Phi=0 [deg];
RP_ref=[0.323293498081869   0.300351794512149]; %C
%}

%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =100;         % lowest number of modes
nMax_u    =100;         % highest number of modes
nMax_step = 1;          % length of the step - modes

N_l    = 200;           % lower number of layers
N_u    = 200;           % upper number of layers
N_step = 1;             % length of the step - layers

tic                   % start time count
RPvec=[];
TPvec=[];

RSvec=[];
TSvec=[];

for nMax =  nMax_l:nMax_step:nMax_u
    
    for N = N_l:N_step:N_u
        tic
        
    %% Relief profile
      d1 = d/N;
      zV = linspace(d1/2,d-d1/2,N);

      t=zeros(2,N);
     xt=zeros(2,N);
     xt_original=zeros(2,N);
     wV=zeros(1,N);
     epst = zeros(2,N);
     
     cosx=zeros(2,N);
     sinx=zeros(2,N);
     x_cs=zeros(1,N);

     epst(1,:)= epsB;
     epst(2,:) = epsW;

    for k=1:N
    t(1,k)=acos(-(2*(zV(k))/d-1))/(2*pi/Lam/cos(Phi));
    t(2,k)=Lam*cos(Phi)-t(1,k);
    
    cosx(1,k)=fcos(t(1,k));
    cosx(2,k)=fcos(t(2,k));
    
    sinx(1,k)=fsin(t(1,k));
    sinx(2,k)=fsin(t(2,k));
    
    xt(1,k)=sec(Phi)*t(1,k)+tan(Phi)*zV(k);
    xt(2,k)=sec(Phi)*t(2,k)+tan(Phi)*zV(k);
    
    xt_original(1,k)=xt(1,k);
    xt_original(2,k)=xt(2,k);
    
    xt(1,k)=xt(1,k)/Lam;
    xt(2,k)=xt(2,k)/Lam;
    if xt(1,k)<1 && xt(2,k)>1
        x0=xt(1,k);
        xt(1,k)=mod(xt(2,k),1);
        xt(2,k)=x0;
        epst(1,k)= epsW;
        epst(2,k) = epsB;
        
        cos0=cosx(1,k);
        cosx(1,k)=cosx(2,k);
        cosx(2,k)=cos0;
       
       sin0=sinx(1,k);
       sinx(1,k)=sinx(2,k);
       sinx(2,k)=sin0;
    
    elseif xt(1,k)>1
        xt(1,k)=mod(xt(1,k),1);
        xt(2,k)=mod(xt(2,k),1);
    end
    end
    %}
     x_cs=(xt(1,:)+xt(2,:))/2;
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(StrucParam.PlotField,'NO')
        [RP,s0V,TP,sSubV] = computeScatMatNVM(lam,thI,epsB,Lam,d,epsS,sinx,cosx,xt,x_cs,epst,nMax,N,StrucParam);
        else
        [RP,s0V,TP,sSubV,H] =computeScatMatNVM_Field(lam,thI,epsB,Lam,d,epsS,sinx,cosx,xt,xt_original,zV,x_cs,epst,nMax,N,StrucParam) ;
        end
 
        % RESULTS
        index_R=find(imag(s0V)==0);
        RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = (abs((RP(index_R)').^2).*s0V(1,index_R))./s0V(1,nMax+1);
        errorP_R((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-RP_ref).^2./RP_ref.^2)));      
        %RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = abs((RS(index_R,nMax+1)').^2).*s0V(1,index_R)./s0V(1,nMax+1);
        %errorS_R((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-RS_ref).^2./RS_ref.^2)));
        
        %{
        index_T=find(imag(sSubV)==0);
        TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  =epsB/epsS* (abs((TP(index_T,nMax+1)').^2).*sSubV(1,index_T))./s0V(1,nMax+1);
        TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = abs((TS(index_T,nMax+1)').^2).*sSubV(1,index_T)./s0V(1,nMax+1);
        errorS_T((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-TS_ref).^2./TS_ref.^2)));
        errorP_T((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-TP_ref).^2./TP_ref.^2)));
        %}
        
              
        c_time((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = toc;
    end
end
       tot_Run_time=sum(sum(c_time));
%% Save results to a file
%filename = 'test.mat';
%save(filename)
%exit
%save lossy_overhang_sinusodial_FMM_FFF_vs_LPEM4.mat
RPvec
RP_ref