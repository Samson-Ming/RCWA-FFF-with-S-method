function [RP,s0V,TP,sSubV,H] =computeScatMatNVM_Field(lam,thI,epsB,Lam,d,epsS,sinx,cosx,xt,xt_original,zV,x_cs,epst,nMax,N,StrucParam) 
%% Inputs are parameters of the grating, incident light and truncation number
%% Outputs are scattering matrices for s- and p-polarization and vector of propagation modes in superstrate


nDim = 2*nMax+1;
k0=2*pi/lam;
q=lam/Lam;
q0=sqrt(epsB)*sin(thI);
I=eye(nDim);
zero=zeros(nDim);
d1 = d/N;
nV=-nMax:nMax;

%% Prepare Factorization Matrices
qV  = q0+nV*q;
kX   = diag(qV);  
    
s0V = sqrt(epsB-qV.^2);
sSubV = sqrt(epsS-qV.^2);

%% W0 and W_Np1 
    W_0_11=I;    
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
    Phi_N=zeros(nDim,nDim,N);
    sVP=zeros(nDim,N);
    sVN=zeros(nDim,N);
  
     %p pol
    %Solve eigenfunction and generate W matrix;
    for iL=1:N   
            epsMn = epsMatrix(epst(1,iL),epst(2,iL),xt(:,iL),nMax);
            etaMn = epsMatrix(1/epst(1,iL),1/epst(2,iL),xt(:,iL),nMax);
            [cosM,sinM] = generateSinCosMat(sinx(:,iL),cosx(:,iL),x_cs(:,iL),nMax); 
                      
            if strcmp(StrucParam.CS,'C_S')
                [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn,StrucParam);
                CP = [-kX*(D\C) I-kX*(D\kX); (A-B*(D\C)) -B*(D\(kX))];
             elseif strcmp(StrucParam.CS,'CC_CS') 
                [A,B,~,D] = generateNVMMat(cosM,sinM,etaMn,epsMn,StrucParam);
                CP = -[kX*D*B,(kX*D*kX-I);(B*D*B-A-etaMn\I),B*D*kX];  
            end
              %CP=double(single(CP));
            [GPH,GNH,GPE,GNE,PP_p,PP_n_inv,sVPos,sVNeg,PP_n] = generatePRPMat(CP,nDim,k0,d1,StrucParam);
            
           W_2_11(:,:,iL)=GPH;
           W_2_12(:,:,iL)=GNH;
           W_2_21(:,:,iL)=GPE;
           W_2_22(:,:,iL)=GNE;
           Phi_P(:,:,iL)=PP_p;
           Phi_N_inv(:,:,iL)=PP_n_inv;
           Phi_N(:,:,iL)=PP_n;
           sVP(:,iL)=sVPos;
           sVN(:,iL)=sVNeg;
    end
%% S matrix: W-t-S
    R_ud=zeros(nDim,nDim,N+1);
    R_du=zeros(nDim,nDim,N+1);
    T_uu=zeros(nDim,nDim,N+1);
    T_dd=zeros(nDim,nDim,N+1);
    
    s_2=zeros(2*nDim,2*nDim,N);
    
    %%{
    s_Np1=([W_Np1_11 -W_2_12(:,:,N) ; W_Np1_21 -W_2_22(:,:,N)])\[W_2_11(:,:,N) -W_Np1_12 ; W_2_21(:,:,N) -W_Np1_22];
   for iL=1:N-1
    s_2(:,:,iL)=([W_2_11(:,:,iL+1) -W_2_12(:,:,iL) ; W_2_21(:,:,iL+1) -W_2_22(:,:,iL)])\[W_2_11(:,:,iL) -W_2_12(:,:,iL+1) ; W_2_21(:,:,iL) -W_2_22(:,:,iL+1)];
   end
   s_0=([W_2_11(:,:,1) -W_0_12 ; W_2_21(:,:,1) -W_0_22])\[W_0_11 -W_2_12(:,:,1) ; W_0_21 -W_2_22(:,:,1)];
   %}
   
   %a_uu(:,:,1)=s_0(1:number_of_orders,1:number_of_orders);
R_ud(:,:,1)=s_0(1:nDim,nDim+1:2*nDim);
%a_du(:,:,1)=s_0(1+number_of_orders:2*number_of_orders,1:number_of_orders);
T_dd(:,:,1)=s_0(nDim+1:2*nDim,nDim+1:2*nDim);

for i=1:1:N
    if i==N
        s_2_l=[I zero;zero Phi_N_inv(:,:,N)]*s_Np1*[Phi_P(:,:,N) zero;zero I];
        
        b_uu=s_2_l(1:nDim,1:nDim);
        b_ud=s_2_l(1:nDim,1+nDim:2*nDim);
        b_du=s_2_l(1+nDim:2*nDim,1:nDim);
        b_dd=s_2_l(nDim+1:2*nDim,nDim+1:2*nDim);
    else
        %s_1_l=[I zero;zero X1(:,:,i)]*s_1(:,:,i)*[X1(:,:,i) zero;zero I];
        s_2_l=diag([I zero;zero Phi_N_inv(:,:,i)])*conj((diag([Phi_P(:,:,i) zero;zero I]))').*s_2(:,:,i);
        
        b_uu=s_2_l(1:nDim,1:nDim);    
        b_ud=s_2_l(1:nDim,1+nDim:2*nDim);    
        b_du=s_2_l(1+nDim:2*nDim,1:nDim);    
        b_dd=s_2_l(nDim+1:2*nDim,nDim+1:2*nDim);
    end
    
    %a_uu(:,:,i+1)=(b_uu/(I-a_ud(:,:,i)*b_du))*a_uu(:,:,i);
    R_ud(:,:,i+1)=b_ud+b_uu*(R_ud(:,:,i)/(I-b_du*R_ud(:,:,i)))*b_dd;
    %a_du(:,:,i+1)=a_du(:,:,i)+a_dd(:,:,i)*(b_du/(I-a_ud(:,:,i)*b_du))*a_uu(:,:,i);
    T_dd(:,:,i+1)=(T_dd(:,:,i)/(I-b_du*R_ud(:,:,i)))*b_dd;    
end

 
    %---------------------------------------------------------
        d_p1=zeros(nDim,1);
        d_p1(nMax+1,1)=1;
        
        u_0=zeros(nDim,1);
%---------------------------------------------------------

        RP=R_ud(:,:,end)*d_p1;
        TP=T_dd(:,:,end)*d_p1; 

%H=0;
%%{
%--------------------------------------------------------- 
d_3=StrucParam.h_upper*d;      %h_upper=StrucParam.h_upper*h_grating;
d_0=StrucParam.h_lower*d;      %h_lower=StrucParam.h_lower*h_grating;
d_tot=d_0+d+ d_3;
a_0=[diag(exp(-1i*k0*sSubV*d_0))*u_0;TP];
a_p1=[RP;diag(exp(-1i*k0*s0V*d_3))*d_p1];  %

a_2=zeros(2*nDim,N);

%%{
a_2(1:nDim,1)=s_0(1:nDim,nDim+1:2*nDim)*(s_0(nDim+1:2*nDim,nDim+1:2*nDim)\TP);

% Pomocí translace volným prostorem provedu transformaci u_1_0 na u_1_1
%u_1_1=X1(:,:,1)*c_m(:,:,N);

% Výpočet c_p(:,:,1) --- z=0

%a_2(1:nDim,N)=s_Np1(1:nDim,1:nDim)\(RP-s_Np1(1:nDim,1+nDim:2*nDim)*d_p1);
u_n_n=s_Np1(1:nDim,1:nDim)\(RP-s_Np1(1:nDim,1+nDim:2*nDim)*d_p1);

a_2(nDim+1:2*nDim,N)=s_Np1(nDim+1:2*nDim,1:nDim)*u_n_n+s_Np1(nDim+1:2*nDim,nDim+1:2*nDim)*d_p1;

d_n_plus_1_n(:,N)=Phi_N_inv(:,:,N)*a_2(nDim+1:2*nDim,N); 
%nepřehledný for-cyklus, nejlépe kouknou do poznámek
if N>=2    
    for i=1:N-1
        a_2(1:nDim,N-i+1)=R_ud(:,:,N-i+1)*d_n_plus_1_n(:,N-i+1); %využívám S-matice --- koeficient R_ud           
        
        inv_s2=inv(s_2(1:nDim,1:nDim,N-i));
        
       u_n_n= inv_s2*(a_2(1:nDim,N-i+1)-s_2(1:nDim,1+nDim:2*nDim,N-i)*a_2(1+nDim:2*nDim,N-i+1));
        
        a_2(nDim+1:2*nDim,N-i)= s_2(1+nDim:2*nDim,1:nDim,N-i)*u_n_n+ s_2(nDim+1:2*nDim,nDim+1:2*nDim,N-i)*d_n_plus_1_n(:,N-i+1);
       
        d_n_plus_1_n(:,N-i)=Phi_N_inv(:,:,N-i)*a_2(nDim+1:2*nDim,N-i); 
    end
end
%}


%Plot field
x=(linspace(0,StrucParam.Number_of_Period*Lam,StrucParam.Resolution_x+1))';

exp_x=(exp(1i*k0*x*qV)).';

% boundary of layers
%{
 boundary=zeros(1,N+1);
 for i=1:1:N
      boundary(i+1)=sum(layer_thickness(1:i))*1e6;    
 end
%}
 boundary=(0:N)*d1;
 
z= (linspace(-d_0,d_tot-d_0,StrucParam.Resolution_z+1))';

d1 = d/N;
boundary=(0:N)*d1;
boundary=[-d_0,boundary,d_tot-d_0];
a_2=[a_0,a_2,a_p1];

H=zeros(StrucParam.Resolution_z+1,StrucParam.Resolution_x+1);

Hymq_u=zeros(nDim,nDim,N+2);
Hymq_u(:,:,1)=I;
Hymq_u(:,:,N+2)=I;
Hymq_u(:,:,2:N+1)=W_2_11(:,:,1:N);

Hymq_d=zeros(nDim,nDim,N+2);
Hymq_d(:,:,1)=I;
Hymq_d(:,:,N+2)=I;
Hymq_d(:,:,2:N+1)=W_2_12(:,:,1:N);
Qu=k0*[sSubV.',sVP,s0V.'];
Qd=k0*[-sSubV.',sVN,-s0V.'];
 
 for z_index=1:StrucParam.Resolution_z+1
        for iL=1:N+2
            if z(z_index)>=boundary(iL) && z(z_index)<=boundary(iL+1)
            Hym=Hymq_u(:,:,iL)*(a_2(1:nDim,iL).*exp(1i*Qu(:,iL)*(z(z_index)-boundary(iL))))+Hymq_d(:,:,iL)*(a_2(nDim+1:2*nDim,iL).*exp(1i*Qd(:,iL)*(z(z_index)-boundary(iL+1))));
            H(z_index,:)=(Hym.')*exp_x;
            end
        end
 end
 
 figure;
 %imagesc(x,z,abs(H));
pcolor(x,z,abs(H));
shading interp
 axis xy
 colorbar
 colormap jet

 hold on
 for ix=1:StrucParam.Number_of_Period
     mid=(xt_original(1,end)+xt_original(2,end))/2;
     xt_original1=[[0;Lam],xt_original,[mid;mid]];
     zV1=[0,zV,d];
     plot(xt_original1-2*Lam+(ix-1)*Lam,zV1,'w','linewidth',1.5);
      plot(xt_original1-Lam+(ix-1)*Lam,zV1,'w','linewidth',1.5);
      plot(xt_original1+(ix-1)*Lam,zV1,'w','linewidth',1.5);
      plot(xt_original1+Lam+(ix-1)*Lam,zV1,'w','linewidth',1.5);
      plot(xt_original1+2*Lam+(ix-1)*Lam,zV1,'w','linewidth',1.5);
 end
 hold off
%}            