%% Inputs are toeplitz matrices of sin and cos functions to profile and 
% toeplitz matrices of permittivities, 
%% Outputs are auxiliary matrices A,B,C,D of NVM
function [A,B,C,D] = generateNVMMat(cosM,sinM,etaMn,epsMn,StrucParam)
     if strcmp(StrucParam.CS,'C_S')
        A=cosM*(etaMn\cosM) + sinM*epsMn*sinM; %Auxiliary factorization matrices
        B=sinM*epsMn*cosM - cosM*(etaMn\sinM);
        C=cosM*epsMn*sinM - sinM*(etaMn\cosM);
        D=sinM*(etaMn\sinM) + cosM*epsMn*cosM;
    elseif strcmp(StrucParam.CS,'CC_CS') 
        I=eye(size(cosM));
        Delta=epsMn-etaMn\I;        
        A=Delta*cosM; %Auxiliary factorization matrices
        B=Delta*sinM;
        D=(epsMn-A)\I;
        C=I;
     end
end



