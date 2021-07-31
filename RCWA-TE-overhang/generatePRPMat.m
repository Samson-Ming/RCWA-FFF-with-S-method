%function [Z_p0,Z_n0,PP_p,PP_n] = generatePRPMat(CP,nDim,k0,d1)
function [GE,GH,PP_p] = generatePRPMat(CP,k0,di,StrucParam)
%% Inputs are CP matrix, truncation number, wavenumber and layer thickness
%% Outputs are admittance and propagation matrices
[GE,sV] = eig(CP);
   sV      = sqrt(diag(sV));
   sVNeg=find(imag(sV)<0);
   sV(sVNeg)=-1*sV(sVNeg);
   GH     = -GE*diag(sV);
   %%
   %%{
    if StrucParam.filtering
    % filtering spurious modes
    sV(real(sV)>StrucParam.threshold)=StrucParam.threshold+1i*1e40;
    end
    %}
    PP_p    = diag(exp(1i*k0*sV*di));  
end