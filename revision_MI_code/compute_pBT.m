function [Syy,f] = compute_pBT(y,maxlag,nfft,fs,lagwindow)
% [Syy,f] = compute_pBT(y,maxlag,nfft,fs,lagwindow)
% perform Blackman-Tukey spectral estimation
% returns power spectrum and corresponding frequencies
% lagwindows available: 
%   't' for a triangular window
% 	'g' for a Gauss window 
% 	'h' for a Hamming window 
% 	'b' for a Blackman window. 

% compute autocorrelation up to maxlag
ryy = xcorr(y-mean(y),maxlag,'biased');
ryy=ryy';
if lagwindow == 't'
   w = triang(2*maxlag+1);
elseif lagwindow == 'g'
   w = gauss(2*maxlag+1);
elseif lagwindow == 'h'
   w = hamming(2*maxlag+1);
elseif lagwindow == 'b'
   w = blackman(2*maxlag+1);   
end
% compute fft
Syy = abs(fft(w.*ryy, nfft));
deltf=fs/nfft;
f=0:deltf:(fs-deltf);