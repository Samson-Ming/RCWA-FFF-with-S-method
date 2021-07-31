function epsM = epsMatrix(epsB, epsW, xt, nMax)
%Toeplitz matrix of permittivity coefficients
%   of a symmetric grating with w=fillfactor

nDim = 2*nMax+1;
kV = 1:2*nMax;

colV = zeros(nDim,1);

colV(1,1) = (1-xt(2)+xt(1))*epsB+(xt(2)-xt(1))*epsW;
%colV(2:(nMax+1),1) = (1i./(2*pi*kV.')).*(exp(-1i*pi*w*kV.')-exp(1i*pi*w*kV.'))*(epsW-epsB);
colV(2:nDim,1) = (1i./(2*pi*kV.')).*(exp(-2i*pi*xt(2)*kV.')-exp(-2i*pi*xt(1)*kV.'))*(epsW-epsB);

rowV = zeros(1,nDim);
rowV(1,1) = (1-xt(2)+xt(1))*epsB+(xt(2)-xt(1))*epsW;
%rowV(1,2:(nMax+1)) = (-1i./(2*pi*kV)).*(exp(1i*pi*w*kV)-exp(-1i*pi*w*kV))*(epsW-epsB);
rowV(1,2:nDim) = (-1i./(2*pi*kV)).*(exp(2i*pi*xt(2)*kV)-exp(2i*pi*xt(1)*kV))*(epsW-epsB);

epsM = toeplitz(colV,rowV);

end
