% planar diffraction by a sinusoidal relif grating
clear all
clc
lam    = 1;
%permW = 1.0;
%permW = 1.72^2;
permW = 1.5^2;
%permW  = (1+5i)^2;   %NOTE that the permitivity is complex!!!
thI    = 0*(pi/180);%+1e-3*(pi/180);
epsB   = 1;
Lam    = 1.5;
d      = 1.5;
epsW   = permW;
epsS   = epsW;
k0     = 2*pi/lam;
q      = lam/Lam;
q0     = sin(thI);

Phi=0*(pi/180);
%Phi=30*(pi/180); 

%% Control parameters
StrucParam.filtering=0; %0: default, no filtering; 1: filtering;; According to Nikolay M. Lyndin, et. al.                                         
StrucParam.threshold=20;  %for eigenvalue filtering; 20: No filtering; 5: used for filtering by experience for highly conducting metal                                      
%% Preallocate fields for reflection amplitudes and errors
errorS_R=zeros(1,1);
errorS_T=zeros(1,1);

%% Reference values for sinusiodal grating, computed by C-Method

%% Reference values for sinusiodal grating, computed by C-Method
%%{
%%n=1.5,theta=0[deg],Phi=14.03 [deg]; C method,
RS_ref=[0.0116110339181890	0.00397254128892009	0.00269456222006578];
TS_ref=[0.0984302860119821	0.476629669959978	0.0994663070054246	0.268814629169694	0.0384145463940903];
%}


%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =13;         % lowest number of modes
nMax_u    =13;         % highest number of modes
nMax_step = 1;          % length of the step - modes

N_l    = 23;           % lower number of layers
N_u    =23;           % upper number of layers
N_step = 2;             % length of the step - layers

tic                   % start time count
%RPvec=[];
%TPvec=[];

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
    %t(1,k)=acos(-(2*(zV(k))/d-1))/(2*pi/Lam/cos(Phi));
    %t(2,k)=Lam*cos(Phi)-t(1,k);
    t(1,k)=3/8*zV(k)+Lam/4;
    t(2,k)=1/8*zV(k)+3/4*Lam;
        
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
    %}
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [RS,s0V,TS,sSubV] = computeScatMatNVM (lam,thI,epsB,Lam,d,epsS,xt,epst,nMax,N,StrucParam);                                                                               
 
        % RESULTS
        index_R=find(imag(s0V)==0);
        %RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = (abs((RP(index_R)').^2).*s0V(1,index_R))./s0V(1,nMax+1);
        %errorP_R((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-RS_ref).^2./RS_ref.^2)));      
        RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = abs((RS(index_R)').^2).*s0V(1,index_R)./s0V(1,nMax+1);
        %errorS_R((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-RS_ref).^2./RS_ref.^2)));
        RS=RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:);
        RS=RS(:);
        %%{
        index_T=find(imag(sSubV)==0);
        TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = abs((TS(index_T)').^2).*sSubV(1,index_T)./s0V(1,nMax+1);
        TS=TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:);
        TS=TS(:);
        %errorS_T((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((TSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-TS_ref).^2./TS_ref.^2)));
        %TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  =epsB/epsS* (abs((TP(index_T)').^2).*sSubV(1,index_T))./s0V(1,nMax+1);       
        %errorP_T((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((TPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-TP_ref).^2./TP_ref.^2)));
        RT=[RS;TS];
        RT_ref=[RS_ref(:);TS_ref(:)];
        %error((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)=abs(RT-RT_ref)./RT_ref;
        error((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)=max(abs(RT-RT_ref)./RT_ref);       
        %}
                      
        c_time((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = toc;
    end
end
%       tot_Run_time=sum(sum(c_time));
%% Save results to a file
       tot_Run_time=sum(sum(c_time))
%       error_study=error(:,:,7);
%% Save results to a file
filename = 'dielectric_overhang_relief_trapezoid_TE_FMM_1percent_N=27_L=23.mat';
%filename = 'dielectric_overhang_relief_trapezoid_TE_FMM_1percent.mat';
save(filename)
%exit
%save lossy_overhang_sinusodial_FMM_FFF_vs_LPEM4.mat
%{
RS
RS_ref
TS
TS_ref
error_study
%}