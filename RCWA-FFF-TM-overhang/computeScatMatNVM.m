function [RP,s0V,TP,sSubV] =computeScatMatNVM(lam,thI,epsB,Lam,d,epsS,sinx,cosx,xt,x_cs,epst,nMax,N,StrucParam)      
%% Inputs are parameters of the grating, incident light and truncation number
%% Outputs are scattering matrices for s- and p-polarization and vector of propagation modes in superstrate
nDim = 2*nMax+1;
k0=2*pi/lam;
q=lam/Lam;
q0=sqrt(epsB)*sin(thI);
I=eye(nDim);
d1 = d/N;

%% Prepare Factorization Matrices
nV=-nMax:nMax;
qV  = q0+nV*q;
kX   = diag(qV);  
    
s0V = sqrt(epsB-qV.^2);
sSubV = sqrt(epsS-qV.^2);
 
%% W0 and W_Np1 
    W_0_11=I;                             % W_1_11=zero,%W_1_21=I
    W_0_12=I;
    W_0_21=diag(sSubV/epsS);
    W_0_22=-W_0_21;
    
    W_Np1_11=I;
    W_Np1_12=I;
    W_Np1_21=diag(s0V/epsB);
    W_Np1_22=-W_Np1_21;
    
    R_ud=zeros(nDim);
    R_du=zeros(nDim);
    T_uu=I;
    T_dd=I;

    %%W2
    W_2_11_i=zeros(nDim,nDim);
    W_2_12_i=zeros(nDim,nDim);
    W_2_21_i=zeros(nDim,nDim);
    W_2_22_i=zeros(nDim,nDim);
    Phi_P_i=zeros(nDim,nDim);
    Phi_N_inv_i=zeros(nDim,nDim);
    
    W_2_11_ip1=zeros(nDim,nDim);
    W_2_12_ip1=zeros(nDim,nDim);
    W_2_21_ip1=zeros(nDim,nDim);
    W_2_22_ip1=zeros(nDim,nDim);
    Phi_P_ip1=zeros(nDim,nDim);
    Phi_N_inv_ip1=zeros(nDim,nDim);
    
%% S matrix: W---S    
    %For the substrate
    epsMn = epsMatrix(epst(1,1),epst(2,1),xt(:,1),nMax);
    etaMn = epsMatrix(1/epst(1,1),1/epst(2,1),xt(:,1),nMax);
   [cosM,sinM] = generateSinCosMat(sinx(:,1),cosx(:,1),x_cs(:,1),nMax); 
   
   if strcmp(StrucParam.CS,'C_S')
   [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn,StrucParam);
   CP = [-kX*(D\C) I-kX*(D\kX); (A-B*(D\C)) -B*(D\(kX))];
   elseif strcmp(StrucParam.CS,'CC_CS') 
   [A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn,StrucParam);
   CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];  
    end
   
   %CP=double(single(CP));
   [GPH,GNH,GPE,GNE,PP_p,PP_n_inv] = generatePRPMat(CP,nDim,k0,d1,StrucParam);
            
   W_2_11_i=GPH;
   W_2_12_i=GNH;
   W_2_21_i=GPE;
   W_2_22_i=GNE;
   Phi_P_i=PP_p;
   Phi_N_inv_i=PP_n_inv;
   
   R_ud_tilde=R_ud;
   T_dd_tilde=T_dd;
   T_uu_tilde=T_uu;
   X1=[W_0_11*T_uu_tilde;W_0_21*T_uu_tilde];
   X2=[-W_2_12_i;-W_2_22_i];
   Z=[W_2_11_i,-W_0_11*R_ud_tilde-W_0_12;...
       W_2_21_i,-W_0_21*R_ud_tilde-W_0_22];
            
   %invZ=inv(Z);         
    matrix1=Z\X1;
    matrix2=Z\X2;
    R_ud=matrix2(1:nDim,:);
    T_dd=T_dd_tilde*matrix2(nDim+1:2*nDim,:);
    T_uu=matrix1(1:nDim,:);
    R_du=R_du+T_dd_tilde*matrix1(nDim+1:2*nDim,:);     
    
    
    %% For layer 2:N
    %Solve eigenfunction and generate W matrix;
    for iL=1:N-1   
            epsMn = epsMatrix(epst(1,iL+1),epst(2,iL+1),xt(:,iL+1),nMax);
            etaMn = epsMatrix(1/epst(1,iL+1),1/epst(2,iL+1),xt(:,iL+1),nMax);
            [cosM,sinM] = generateSinCosMat(sinx(:,iL+1),cosx(:,iL+1),x_cs(:,iL+1),nMax); 
            
             if strcmp(StrucParam.CS,'C_S')
                [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn,StrucParam);
                CP = [-kX*(D\C) I-kX*(D\kX); (A-B*(D\C)) -B*(D\(kX))];
             elseif strcmp(StrucParam.CS,'CC_CS') 
                [A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn,StrucParam);
                CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];  
              end
            
            
            %CP=double(single(CP));
            [GPH,GNH,GPE,GNE,PP_p,PP_n_inv] = generatePRPMat(CP,nDim,k0,d1,StrucParam);
            
           W_2_11_ip1=GPH;
           W_2_12_ip1=GNH;
          W_2_21_ip1=GPE;
          W_2_22_ip1=GNE;
          Phi_P_ip1=PP_p;
          Phi_N_inv_ip1=PP_n_inv;
          
          R_ud_tilde=diag(Phi_P_i)*((diag(Phi_N_inv_i)).').*R_ud;   %faster; in fact: Phi_p*R_ud*inv(Phi_n);
          T_dd_tilde=T_dd*sparse(Phi_N_inv_i);
           T_uu_tilde=sparse(Phi_P_i)*T_uu;
           X1=[W_2_11_i*T_uu_tilde;W_2_21_i*T_uu_tilde];
           X2=[-W_2_12_ip1;-W_2_22_ip1];
           
           Z=[W_2_11_ip1,-W_2_11_i*R_ud_tilde-W_2_12_i;...
                W_2_21_ip1,-W_2_21_i*R_ud_tilde-W_2_22_i];
           
        %invZ=inv(Z);         
        matrix1=Z\X1;
        matrix2=Z\X2;
            
        R_ud=matrix2(1:nDim,:);
        T_dd=T_dd_tilde*matrix2(nDim+1:2*nDim,:);
        T_uu=matrix1(1:nDim,:);
        R_du=R_du+T_dd_tilde*matrix1(nDim+1:2*nDim,:);           
        
        W_2_11_i=W_2_11_ip1;
        W_2_12_i=W_2_12_ip1;
        W_2_21_i=W_2_21_ip1;
       W_2_22_i=W_2_22_ip1;
       Phi_P_i=Phi_P_ip1;
       Phi_N_inv_i=Phi_N_inv_ip1;
    end
    
    % For layer N
   R_ud_tilde=diag(Phi_P_i)*((diag(Phi_N_inv_i)).').*R_ud;   %faster; in fact: Phi_p*R_ud*inv(Phi_n);
   T_dd_tilde=T_dd*sparse(Phi_N_inv_i);
   T_uu_tilde=sparse(Phi_P_i)*T_uu;
   X1=[W_2_11_i*T_uu_tilde;W_2_21_i*T_uu_tilde];
   X2=[-W_Np1_12;-W_Np1_22];
   Z=[W_Np1_11,-W_2_11_i*R_ud_tilde-W_2_12_i;...
        W_Np1_21,-W_2_21_i*R_ud_tilde-W_2_22_i];
         
    %invZ=inv(Z);         
        matrix1=Z\X1;
        matrix2=Z\X2;
            
        R_ud=matrix2(1:nDim,:);
        T_dd=T_dd_tilde*matrix2(nDim+1:2*nDim,:);
        T_uu=matrix1(1:nDim,:);
        R_du=R_du+T_dd_tilde*matrix1(nDim+1:2*nDim,:); 
        
    % R and T
 %---------------------------------------------------------
        d_p1=zeros(nDim,1);
        d_p1(nMax+1,1)=1;
%---------------------------------------------------------

        RP=R_ud(:,:,end)*d_p1;
        TP=T_dd(:,:,end)*d_p1;  
        