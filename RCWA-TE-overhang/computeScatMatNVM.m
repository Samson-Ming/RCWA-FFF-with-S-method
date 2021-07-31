function [RP,s0V,TP,sSubV] =computeScatMatNVM (lam,thI,epsB,Lam,d,epsS,xt,epst,nMax,N,StrucParam)      
%% Inputs are parameters of the grating, incident light and truncation number
%% Outputs are scattering matrices for s- and p-polarization and vector of propagation modes in superstrate
nDim = 2*nMax+1;
k0=2*pi/lam;
q=lam/Lam;
q0=sin(thI);
I=eye(nDim);
d1 = d/N;

%% Prepare Factorization Matrices

nV=-nMax:nMax;
qV  = sqrt(epsB)*q0+nV*q;
kX   = diag(qV);  
    
s0V = sqrt(epsB-qV.^2);
sSubV = sqrt(epsS-qV.^2);
 
%% W0 and W_Np1 
    W_0_1=I;                             % W_1_11=zero,%W_1_21=I
    W_0_2=-diag(sSubV);
    W_Np1_1=I;
    W_Np1_2=-diag(s0V);
    
    R_ud=zeros(nDim);
    R_du=zeros(nDim);
    T_uu=I;
    T_dd=I;

    %%W2
    W_2_1_i=zeros(nDim,nDim);
    W_2_2_i=zeros(nDim,nDim);
    Phi_P_i=zeros(nDim,nDim);
    
    W_2_1_ip1=zeros(nDim,nDim);
    W_2_2_ip1=zeros(nDim,nDim);
    Phi_P_ip1=zeros(nDim,nDim);
    
%% S matrix: W---S    
    %For the substrate
    epsMn = epsMatrix(epst(1,1),epst(2,1),xt(:,1),nMax);
   
   CP=epsMn-kX^2;

   [GE,GH,PP_p] = generatePRPMat(CP,k0,d1,StrucParam);
            
   W_2_1_i=GE;
   W_2_2_i=GH;
   Phi_P_i=PP_p;
   
   R_ud_tilde=R_ud;
   T_dd_tilde=T_dd;
   T_uu_tilde=T_uu;
   
   Q1=W_2_1_i\W_0_1;
   Q2=W_2_2_i\W_0_2;
   F=Q1*(I+R_ud_tilde);
   G=Q2*(I-R_ud_tilde);
   
   tau=(F+G)\I;
   
    R_ud=I-2*G*tau;
    T_dd=2*T_dd_tilde*tau;
    T_uu=(F*tau*Q2+G*tau*Q1)*T_uu_tilde;
    R_du=R_du+T_dd_tilde*tau*(Q2-Q1)*T_uu_tilde;
    
    
    %% For layer 1:N-1
    %Solve eigenfunction and generate W matrix;
    for iL=1:N-1   
            epsMn = epsMatrix(epst(1,iL+1),epst(2,iL+1),xt(:,iL+1),nMax);
            
           CP=epsMn-kX^2;

          [GE,GH,PP_p] = generatePRPMat(CP,k0,d1,StrucParam);
            
           W_2_1_ip1=GE;
          W_2_2_ip1=GH;
          Phi_P_ip1=PP_p;
          
          R_ud_tilde=diag(Phi_P_i)*((diag(Phi_P_i)).').*R_ud;   %faster; in fact: Phi_p*R_ud*inv(Phi_n);
          T_dd_tilde=T_dd*sparse(Phi_P_i);
          T_uu_tilde=sparse(Phi_P_i)*T_uu;
          
         Q1=W_2_1_ip1\W_2_1_i;
         Q2=W_2_2_ip1\W_2_2_i;
         F=Q1*(I+R_ud_tilde);
         G=Q2*(I-R_ud_tilde);
   
        tau=(F+G)\I;
   
        R_ud=I-2*G*tau;
        T_dd=2*T_dd_tilde*tau;
        T_uu=(F*tau*Q2+G*tau*Q1)*T_uu_tilde;
        R_du=R_du+T_dd_tilde*tau*(Q2-Q1)*T_uu_tilde;
                          
        
        W_2_1_i=W_2_1_ip1;
        W_2_2_i=W_2_2_ip1;
       Phi_P_i=Phi_P_ip1;
    end
    
    % For layer N
   R_ud_tilde=diag(Phi_P_i)*((diag(Phi_P_i)).').*R_ud;   %faster; in fact: Phi_p*R_ud*inv(Phi_n);
   T_dd_tilde=T_dd*sparse(Phi_P_i);
   T_uu_tilde=sparse(Phi_P_i)*T_uu;
   
   Q1=W_Np1_1\W_2_1_i;
   Q2=W_Np1_2\W_2_2_i;
   F=Q1*(I+R_ud_tilde);
   G=Q2*(I-R_ud_tilde);
   
    tau=(F+G)\I;
   
   R_ud=I-2*G*tau;
   T_dd=2*T_dd_tilde*tau;
   T_uu=(F*tau*Q2+G*tau*Q1)*T_uu_tilde;
   R_du=R_du+T_dd_tilde*tau*(Q2-Q1)*T_uu_tilde;
        
    % R and T
 %---------------------------------------------------------
        d_n1=zeros(nDim,1);
        d_n1(nMax+1,1)=1;
%---------------------------------------------------------

        RP=R_ud(:,:,end)*d_n1;
        TP=T_dd(:,:,end)*d_n1;  
        