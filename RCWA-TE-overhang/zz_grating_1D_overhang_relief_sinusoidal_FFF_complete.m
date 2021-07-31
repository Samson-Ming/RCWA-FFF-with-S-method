% planar diffraction by a sinusoidal relif grating
clear all
clc
lam    = 1;
%permW = 1.0;
%permW = 1.72^2;
%permW = 3^2;
permW  = (1+5i)^2;   %NOTE that the permitivity is complex!!!
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
StrucParam.filtering=0; %0: default, no filtering; 1: filtering;; According to Nikolay M. Lyndin, et. al.                                         
StrucParam.threshold=20;  %for eigenvalue filtering; 20: No filtering; 5: used for filtering by experience for highly conducting metal                                      
%% Preallocate fields for reflection amplitudes and errors
errorS_R=zeros(1,1);
errorS_T=zeros(1,1);

%% Reference values for sinusiodal grating, computed by C-Method

%{
%%n=1.72,theta=0 [deg], Phi=30 [deg];
RS_ref=[0.010174769316324];
TS_ref=[0.726658666521634	0.229551910187410	0.033614676815176];
%}

%{
%%n=1.72,theta=15[deg],Phi=30 [deg];
RS_ref=[0.020351069199529   0.002595456005422];
TS_ref=[0.465951079541316   0.429741994958952   0.081360436389336];
%}

%{
%%n=1.72,theta=15[deg],Phi=30 [deg]; LPEM
RS_ref=[0.033479360118408   0.003225124752801];
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
RS_ref=[0.271923709974736   0.502654339036197];
%FEM
RS_ref=[0.27262164737017386	0.5021319693110445];
%Stardard FFM:RCWA-1D+S method:N=401,M=200;
RS_ref=[0.271998824578871 0.502563618626655];
%Stardard FFM:Reticolo+modified T method:N=401,M=200;
RS_ref=[0.271998824578443 0.502563618627200];

%FFM-FFF:N=401,M=200;
RS_ref=[0.271998824578812 0.502563618626977];  %Staircase_CS

%%FFM-FFF-new:N=101,M=300;
RS_ref=[0.271969795841654 0.502621639352014];  %Staircase

%FFM-FFF-new:N=61,M=200;
RS_ref=[0.271982260224289 0.502530746717528];  %Staircase

%C method
RS_ref=[0.271937696305567   0.502697416029626];
%}

%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =200;         % lowest number of modes
nMax_u    =200;         % highest number of modes
nMax_step = 1;          % length of the step - modes

N_l    = 200;           % lower number of layers
N_u    =200;           % upper number of layers
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
     

     epst(1,:)= epsB;
     epst(2,:) = epsW;

    for k=1:N
    t(1,k)=acos(-(2*(zV(k))/d-1))/(2*pi/Lam/cos(Phi));
    t(2,k)=Lam*cos(Phi)-t(1,k);
        
    xt(1,k)=sec(Phi)*t(1,k)+tan(Phi)*zV(k);
    xt(2,k)=sec(Phi)*t(2,k)+tan(Phi)*zV(k);
       
    xt(1,k)=xt(1,k)/Lam;
    xt(2,k)=xt(2,k)/Lam;
    if xt(1,k)<1 && xt(2,k)>1
        x0=xt(1,k);
        xt(1,k)=mod(xt(2,k),1);
        xt(2,k)=x0;
        epst(1,k)= epsW;
        epst(2,k) = epsB;
       
    elseif xt(1,k)>1
        xt(1,k)=mod(xt(1,k),1);
        xt(2,k)=mod(xt(2,k),1);
    end
    end
end
    %}
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [RS,s0V,TS,sSubV] = computeScatMatNVM (lam,thI,epsB,Lam,d,epsS,xt,epst,nMax,N,StrucParam);                                                                               
 
        % RESULTS
        index_R=find(imag(s0V)==0);
        %RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = (abs((RP(index_R)').^2).*s0V(1,index_R))./s0V(1,nMax+1);
        %errorP_R((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-RS_ref).^2./RS_ref.^2)));      
        RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = abs((RS(index_R)').^2).*s0V(1,index_R)./s0V(1,nMax+1);
        errorS_R((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-RS_ref).^2./RS_ref.^2)));
        
        %{
        index_T=find(imag(sSubV)==0);
        TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  =epsB/epsS* (abs((TP(index_T)').^2).*sSubV(1,index_T))./s0V(1,nMax+1);
        TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = abs((TS(index_T)').^2).*sSubV(1,index_T)./s0V(1,nMax+1);
        errorS_T((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-TS_ref).^2./TS_ref.^2)));
        errorP_T((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-TP_ref).^2./TP_ref.^2)));
        %}
        
              
        c_time((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = toc;
    end
%end
       tot_Run_time=sum(sum(c_time));
%% Save results to a file
%filename = 'test.mat';
%save(filename)
%exit
%save lossy_overhang_sinusodial_FMM_FFF_vs_LPEM4.mat
RSvec