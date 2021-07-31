function [RP,s0V,TP,sSubV] = computeScatMatNVM_Field (lam,thI,epsB,Lam,d,epsW,epsS,Phi,fsin,fcos,nMax,N)
%% Inputs are parameters of the grating, incident light and truncation number
%% Outputs are scattering matrices for s- and p-polarization and vector of propagation modes in superstrate


nDim = 2*nMax+1;
k0=2*pi/lam;
q=lam/Lam;
q0=sin(thI);
I=eye(nDim);


%% Relief profile

d1 = d/N;
zV = linspace(d1/2,d-d1/2,N);

xt=zeros(2,N);
xt_original=zeros(2,N);
wV=zeros(1,N);
epst = zeros(2,N);

epst(1,:)= epsB;
epst(2,:) = epsW;

%% NOTE that z is downward, opposite to usall??????
for k=1:N
    xt(1,k)=acos(-(2*(zV(k))/d-1))/(2*pi/Lam/cos(Phi));
    xt(2,k)=Lam*cos(Phi)-xt(1,k);
    
    xt(1,k)=sec(Phi)*xt(1,k)+tan(Phi)*zV(k);
    xt(2,k)=sec(Phi)*xt(2,k)+tan(Phi)*zV(k);
end
%xt=fliplr(xt);
xt=xt/Lam;

%wV = (1/pi)*acos(1-2*zV/d); % vector of fill factors ... sinusoidal profile
%wV = 0.5*ones(1,N);


nV=-nMax:nMax;

%% Prepare Factorization Matrices
m=14; %orders of FFT;
[cosM,sinM] = generateSinCosMat(fcos,fsin, Lam,Phi, nDim,m);
  
qV  = sqrt(epsB)*q0+nV*q;
kX   = diag(qV);  
    
s0V = sqrt(epsB-qV.^2);
sSubV = sqrt(epsS-qV.^2);
 
%% FFF Superstrate medium
%{ 
epsMn = epsMatrix(epsB,epsB,[0;0],nMax);  % "n"="next"
    etaMn = epsMatrix(1/epsB,1/epsB,[0;0],nMax);
    kX   = diag(qV);  
%% Assemble the main matrix, p-polarization
    [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
    CP0 = [-kX0*(D\C) I-kX0*(D\kX0); (A-B*(D\C)) -B*(D\(kX0))];   %NO SENSE???
%% Generate propagation matrices in superstrate, p-polarization
 %   [Z_p0,Z_n0,~,~] = generatePRPMat(CP0,nDim,k0,0);
%}
%% W0 and W_Np1 
    W_0_11=I;                             % W_1_11=zero,%W_1_21=I
    W_0_12=I;
    W_0_21=diag(sSubV/epsS);
    W_0_22=-W_0_21;
    
    W_Np1_11=I;
    W_Np1_12=I;
    W_Np1_21=diag(s0V/epsB);
    W_Np1_22=-W_Np1_21;
    
    W_2_11=zeros(nDim,nDim,N);
    W_2_12=zeros(nDim,nDim,N);
    W_2_21=zeros(nDim,nDim,N);
    W_2_22=zeros(nDim,nDim,N);
    Phi_P=zeros(nDim,nDim,N);
    Phi_N_inv=zeros(nDim,nDim,N);
    
     %p pol
    %Solve eigenfunction and generate W matrix;
    for iL=1:N   
            epsMn = epsMatrix(epst(1,iL),epst(2,iL),xt(:,iL),nMax);
            etaMn = epsMatrix(1/epst(1,iL),1/epst(2,iL),xt(:,iL),nMax);
                      
            [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn);
            CP = [-kX*(D\C) I-kX*(D\kX); (A-B*(D\C)) -B*(D\(kX))];          
            [GPH,GNH,GPE,GNE,PP_p,PP_n_inv] = generatePRPMat(CP,nDim,k0,d1);
            
           W_2_11(:,:,iL)=GPH;
           W_2_12(:,:,iL)=GNH;
           W_2_21(:,:,iL)=GPE;
           W_2_22(:,:,iL)=GNE;
           Phi_P(:,:,iL)=PP_p;
           Phi_N_inv(:,:,iL)=PP_n_inv;
    end
%% S matrix: W---S
    R_ud=zeros(nDim);
    R_du=zeros(nDim);
    T_uu=I;
    T_dd=I;
    
    for iL=0:N
        if iL==0
            R_ud_tilde=R_ud;
            T_dd_tilde=T_dd;
            T_uu_tilde=T_uu;
            X1=[W_0_11*T_uu_tilde;W_0_21*T_uu_tilde];
            X2=[-W_2_12(:,:,iL+1);-W_2_22(:,:,iL+1)];
            Z=[W_2_11(:,:,iL+1),-W_0_11*R_ud_tilde-W_0_12;...
                W_2_21(:,:,iL+1),-W_0_21*R_ud_tilde-W_0_22];
            
            %invZ=inv(Z);         
            matrix1=Z\X1;
            matrix2=Z\X2;
            R_ud=matrix2(1:nDim,:);
            T_dd=T_dd_tilde*matrix2(nDim+1:2*nDim,:);
            T_uu=matrix1(1:nDim,:);
            R_du=R_du+T_dd_tilde*matrix1(nDim+1:2*nDim,:);                  
             
        elseif iL==N          
        R_ud_tilde=diag(Phi_P(:,:,iL))*((diag(Phi_N_inv(:,:,iL))).').*R_ud;   %faster; in fact: Phi_p*R_ud*inv(Phi_n);
        T_dd_tilde=T_dd*sparse(Phi_N_inv(:,:,iL));
        T_uu_tilde=sparse(Phi_P(:,:,iL))*T_uu;
        X1=[W_2_11(:,:,iL)*T_uu_tilde;W_2_21(:,:,iL)*T_uu_tilde];
        X2=[-W_Np1_12;-W_Np1_22];
        Z=[W_Np1_11,-W_2_11(:,:,iL)*R_ud_tilde-W_2_12(:,:,iL);...
             W_Np1_21,-W_2_21(:,:,iL)*R_ud_tilde-W_2_22(:,:,iL)];
         
        %invZ=inv(Z);         
        matrix1=Z\X1;
        matrix2=Z\X2;
            
        R_ud=matrix2(1:nDim,:);
        T_dd=T_dd_tilde*matrix2(nDim+1:2*nDim,:);
        T_uu=matrix1(1:nDim,:);
        R_du=R_du+T_dd_tilde*matrix1(nDim+1:2*nDim,:); 
        
        else
            R_ud_tilde=diag(Phi_P(:,:,iL))*((diag(Phi_N_inv(:,:,iL))).').*R_ud;   %faster; in fact: Phi_p*R_ud*inv(Phi_n);
            T_dd_tilde=T_dd*sparse(Phi_N_inv(:,:,iL));
           T_uu_tilde=sparse(Phi_P(:,:,iL))*T_uu;
           X1=[W_2_11(:,:,iL)*T_uu_tilde;W_2_21(:,:,iL)*T_uu_tilde];
           X2=[-W_2_12(:,:,iL+1);-W_2_22(:,:,iL+1)];
           
           Z=[W_2_11(:,:,iL+1),-W_2_11(:,:,iL)*R_ud_tilde-W_2_12(:,:,iL);...
                W_2_21(:,:,iL+1),-W_2_21(:,:,iL)*R_ud_tilde-W_2_22(:,:,iL)];
           
        %invZ=inv(Z);         
        matrix1=Z\X1;
        matrix2=Z\X2;
            
        R_ud=matrix2(1:nDim,:);
        T_dd=T_dd_tilde*matrix2(nDim+1:2*nDim,:);
        T_uu=matrix1(1:nDim,:);
        R_du=R_du+T_dd_tilde*matrix1(nDim+1:2*nDim,:);               
        end
    end
 
    %---------------------------------------------------------
        d_n1=zeros(nDim,1);
        d_n1(nMax+1,1)=1;
%---------------------------------------------------------

        RP=R_ud(:,:,end)*d_n1;
        TP=T_dd(:,:,end)*d_n1;  
    
    