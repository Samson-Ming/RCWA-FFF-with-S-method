% planar diffraction by a sinusoidal relif grating
clear all
clc
lam    = 1;
%permW = 1.0;
%permW = 1.72^2;
%permW = 3^2;
permW  = (1+5i)^2;   %NOTE that the permitivity is complex!!!
thI    = 45*(pi/180);%+1e-3*(pi/180);
epsB   = 1;
Lam    = 1;
d      = 0.5;
epsW   = permW;
epsS   = epsW;
k0     = 2*pi/lam;
q      = lam/Lam;
q0     = sin(thI);

%Phi=0*(pi/180);
Phi=60*(pi/180); 

%% Control parameters
StrucParam.filtering=0; %0: default, no filtering; 1: filtering;; According to Nikolay M. Lyndin, et. al.                                         
StrucParam.threshold=20;  %for eigenvalue filtering; 20: No filtering; 5: used for filtering by experience for highly conducting metal                                      
%% Preallocate fields for reflection amplitudes and errors
errorS_R=zeros(1,1);
errorS_T=zeros(1,1);

%% Reference values for sinusiodal grating, computed by C-Method


%C method
RS_ref=[0.344941437878085,0.534255315313941];
%}

%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =5;         % lowest number of modes
nMax_u    =5;         % highest number of modes
nMax_step = 1;          % length of the step - modes

N_l    = 15;           % lower number of layers
N_u    =15;           % upper number of layers
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
   x1=mod(xt(1,k),1);
     x2=mod(xt(2,k),1);
     if x1>x2
         xt(1,k)=x2;
         xt(2,k)=x1;
         
         epst(1,k)= epsW;
        epst(2,k) = epsB;
                
     else
         xt(1,k)=x1;
         xt(2,k)=x2;
     end
    end

    %}
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [RS,s0V,TS,sSubV] = computeScatMatNVM (lam,thI,epsB,Lam,d,epsS,xt,epst,nMax,N,StrucParam);                                                                               
 
        % RESULTS
        index_R=find(imag(s0V)==0);
        %RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)  = (abs((RP(index_R)').^2).*s0V(1,index_R))./s0V(1,nMax+1);
        %errorP_R((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) =1/numel(RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:))*sum(sum(sum((RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)-RS_ref).^2./RS_ref.^2)));      
        RS= abs((RS(index_R)').^2).*s0V(1,index_R)./s0V(1,nMax+1);
       RSvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:)=RS;
       RS=RS(:); 
              
        RT=RS;
        RT_ref=RS_ref(:);
        %error_mean((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)=mean(abs(RT-RT_ref)./RT_ref);
        %error_max((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)=max(abs(RT-RT_ref)./RT_ref);
        error((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)=max(abs(RT-RT_ref)./RT_ref);
              
        c_time((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = toc;
    end
end
       tot_Run_time=sum(sum(c_time))
%% Save results to a file
%% Save results to a file
%filename = 'metal_overhang_relief_sinusodial_TE_1percent.mat';
filename = 'metal_overhang_relief_sinusodial_h=0.5lam_Phi=60_TE_1percent_N=11_L=15.mat';
%filename = 'Near_fileld_metal_overhang_relief_sinusodial_h=0.5lam_Phi=60_TE_1percent_N=201_L=300.mat';
save(filename)
%RS