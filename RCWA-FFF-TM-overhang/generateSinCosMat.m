function [cosM,sinM] = generateSinCosMat(sinx,cosx,x_cs,nMax)

%% Inputs are functions of tangent and normal to profile, outputs are Topelitz matrices generated from its fourier coefficients 
%{
nDim = 2*nMax+1;
kV = 1:2*nMax;

cosM = zeros(1,nDim);
cosM(1,1) =x_cs*cosx(1)+(1-x_cs)*cosx(2);
cosM(1,2:nDim) = (-1i./(2*pi*kV)).*(1-exp(2i*pi*x_cs*kV))*(cosx(2)-cosx(1));

sinM = zeros(1,nDim);
sinM(1,1) =x_cs*sinx(1)+(1-x_cs)*sinx(2);
sinM(1,2:nDim) = (-1i./(2*pi*kV)).*(1-exp(2i*pi*x_cs*kV))*(sinx(2)-sinx(1));

cosM = toeplitz(cosM);
sinM = toeplitz(sinM);
%}

if size(x_cs,1)==1
    x_cs=[x_cs*0,x_cs];
end
nDim = 2*nMax+1;
kV = 1:2*nMax;

cosM = zeros(1,nDim);
cosM(1,1) =(x_cs(2)-x_cs(1))*cosx(1)+(1-(x_cs(2)-x_cs(1)))*cosx(2);
cosM(1,2:nDim) = (-1i./(2*pi*kV)).*(exp(2i*pi*x_cs(1)*kV)-exp(2i*pi*x_cs(2)*kV))*(cosx(2)-cosx(1));

sinM = zeros(1,nDim);
sinM(1,1) =(x_cs(2)-x_cs(1))*sinx(1)+(1-(x_cs(2)-x_cs(1)))*sinx(2);
sinM(1,2:nDim) = (-1i./(2*pi*kV)).*(exp(2i*pi*x_cs(1)*kV)-exp(2i*pi*x_cs(2)*kV))*(sinx(2)-sinx(1));

cosM = toeplitz(cosM);
sinM = toeplitz(sinM);

end

