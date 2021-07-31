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
StrucParam.CS='CC_CS'; 
                                      %'CC_CS': default,cos^2 and cos*sin: according to Popov, et.al. , and Habib Mohamad et. al.
                                      %'C_S': cos and sin: according to T.V.
                                                                          
StrucParam.filtering=0; %0: default, no filtering; 1: filtering; According to Nikolay M. Lyndin, et. al.                     
StrucParam.threshold=20;  %for eigenvalue
%Near field plot;
StrucParam.PlotField='YES';  %'NO':default,not plot field;
                                        %'YES',plot near field, quite slow.
StrucParam.h_upper=0.1;      %h_upper=StrucParam.h_upper*h_grating;
StrucParam.h_lower=0.1;      %h_lower=StrucParam.h_lower*h_grating;
StrucParam.Number_of_Period=1;
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
%C method
RP_ref=[0.753199898155597,0.0400481471266204];
%}

%% Set truncation parameters
c_time    = zeros(1,1); % computation time

nMax_l    =100;         % lowest number of modes
nMax_u    =100;         % highest number of modes
nMax_step =1;          % length of the step - modes

N_l    =300;           % lower number of layers, 400
N_u    =300;          % upper number of layers
N_step = 1;             % length of the step -layers

tic                   % start time count
RPvec=[];
TPvec=[];

RSvec=[];
TSvec=[];

for nMax =  nMax_l:nMax_step:nMax_u
    
    for N = N_l:N_step:N_u
        tic
        nMax
        N
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
     %x_cs=zeros(1,N);

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
    %{
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
    %}
     %%{
     x1=mod(xt(1,k),1);
     x2=mod(xt(2,k),1);
     if x1>x2
         xt(1,k)=x2;
         xt(2,k)=x1;
         
         epst(1,k)= epsW;
        epst(2,k) = epsB;
        
        cos0=cosx(1,k);
        cosx(1,k)=cosx(2,k);
        cosx(2,k)=cos0;
       
       sin0=sinx(1,k);
       sinx(1,k)=sinx(2,k);
       sinx(2,k)=sin0;
     else
         xt(1,k)=x1;
         xt(2,k)=x2;
     end
    %}
    end

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
        RP=RPvec((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step,:);
        RP=RP(:); 
        
        RT=RP;
        RT_ref=RP_ref(:);
        %error_mean((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)=mean(abs(RT-RT_ref)./RT_ref);
        %error_max((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)=max(abs(RT-RT_ref)./RT_ref);
        error((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step)=max(abs(RT-RT_ref)./RT_ref);
              
        c_time((nMax-nMax_l+nMax_step)/nMax_step,(N-N_l+N_step)/N_step) = toc;
    end
end
       tot_Run_time=sum(sum(c_time))

%% Save results to a file
%filename = 'metal_overhang_relief_sinusodial_TM_1percent.mat';
%filename = 'metal_overhang_relief_sinusodial_h=0.5lam_Phi=60_TM_1percent_N=201_L=300.mat';
filename = 'Near_fileld_metal_overhang_relief_sinusodial_h=0.5lam_Phi=60_TM_1percent_N=201_L=300.mat';
save(filename)