%function [Z_p0,Z_n0,PP_p,PP_n] = generatePRPMat(CP,nDim,k0,d1)
function [GPH,GNH,GPE,GNE,PP_p,PP_n_inv,sVPos,sVNeg,PP_n] = generatePRPMat(CP,nDim,k0,di,StrucParam)
%% Inputs are CP matrix, truncation number, wavenumber and layer thickness
%% Outputs are admittance and propagation matrices

    [GP,sV] = eig(CP);
    sV      = diag(sV);
    
    %% Partion
    %{
    sVPsort = ((real(sV)+imag(sV))>0);
     sVNsort = ~sVPsort;
     %}
    
%%{
    sVPsort = find((imag(sV))>0);
   
    sVNsort =find((imag(sV))<0);
    Index=(1:2*nDim)';
    sVRsort=setdiff(Index,[sVPsort ;sVNsort]);
    
    if length(sVPsort)>nDim
        [~,indexP]=sort(imag(sV(sVPsort)),'descend');
        sVNsort=[sVNsort;sVRsort;sVPsort(indexP(nDim+1:end))];
        sVPsort=sVPsort(indexP(1:nDim));
    end
    
    if length(sVNsort)>nDim
        [~,indexN]=sort(imag(sV(sVNsort)),'ascend');
        sVPsort=[sVPsort;sVRsort;sVNsort(indexN(nDim+1:end))];
        sVNsort=sVNsort(indexN(1:nDim));
    end
        
   
    
    if (length(sVPsort)<=nDim) && (length(sVNsort)<=nDim) 
        if length(sVPsort)<nDim
            NN=length(sVPsort);
            sVPsort=[sVPsort;sVRsort(1:(nDim-NN))];   
            sVNsort=setdiff(Index,sVPsort);
        elseif length(sVNsort)<nDim
            sVNsort=[sVNsort;sVRsort];
        end
    end
    %}
    
   %% Assembling;
    sVPos   = sV(sVPsort);
    sVNeg   = sV(sVNsort);
    
    %%{
    if StrucParam.filtering
    % filtering spurious modes
    sVPos(real(sVPos)>StrucParam.threshold)=StrucParam.threshold+1i*1e40;
    sVNeg(real(sVNeg)<-StrucParam.threshold)=-StrucParam.threshold-1i*1e40;
    end
    %}
    
    GPE     = GP(1:nDim,sVPsort);
    GPH     = GP(nDim+1:2*nDim,sVPsort);
    GNE     = GP(1:nDim,sVNsort);  
    GNH     = GP(nDim+1:2*nDim,sVNsort);
    %Z_p0    = GPE/GPH;
    %Z_n0    = GNE/GNH;
    PP_p    = diag(exp(1i*k0*sVPos*di));
    PP_n    = diag(exp(1i*k0*sVNeg*di)); 
    PP_n_inv    = diag(exp(-1i*k0*sVNeg*di));    
end