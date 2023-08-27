function [ prate, pcount, pratetime ] = neural_filtering( filter,stimul,alpha,beta,sr,countingwindowlim )

% get filter and stimulus
S=stimul;
F=filter;

% convolve to compute filter response ad apply nonlinearlity
fout=alpha.*((max(S'*F,zeros(size(S,2),1))).^beta)';

% pad filter response at the beginning
prate=[fout(1)*ones(10,1);fout'];

% get filter response time vector
pratetime=(1/sr).*(1:size(S,2))-(1/sr)*25;

% integrate over counting window
countingwindowsamples=and(pratetime>countingwindowlim(1),pratetime<countingwindowlim(2));
pcount=nansum(prate(countingwindowsamples));

end

