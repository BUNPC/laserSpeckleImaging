function [fftPow,fftPhase,f]=getFFT(data,period,fftN,procType)
data=single(data);
if strcmp(procType,'gpu')
    data=gpuArray(data);
end
data=data-mean(data,3);
f = 1./period.*(0:(fftN/2))/fftN;
fftRes = fft(data,fftN,3);
fftRes = fftRes(:,:,1:fftN/2+1);
fftPow = abs(fftRes./fftN/2);
fftPow(:,:,2:end-1)= 2.*fftPow(:,:,2:end-1);
fftPhase=angle(fftRes);
if strcmp(procType,'gpu')
    fftPow=gather(fftPow);
    fftPhase=gather(fftPhase);
end
end
