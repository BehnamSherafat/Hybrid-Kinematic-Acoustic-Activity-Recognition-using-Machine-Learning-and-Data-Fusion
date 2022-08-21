function [A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S, Start, End, fs, FFTSIZE, noverlap, WINDOWSIZE)
%        1-513, 514,515,516, 517,   518,    ,519, 520
%% Extracting 513 Decibel Features
s1 = S(ceil(Start*fs+1):ceil(End*fs));
S_1 = spectrogram(s1, hann(WINDOWSIZE), noverlap, FFTSIZE, 'yaxis');
A_Spec = abs(S_1'); % Normalization
A_SpecM = mean(A_Spec, 2);

%% Calculating RMS and Short time Energy for each window
nofw = floor((length(s1)-noverlap)/noverlap);
A_RMS =zeros(nofw,1);
A_E = zeros(nofw,1);
for m = 1:nofw
    y = s1(round((m-1)*noverlap+1):min(length(s1), floor((m-1)*noverlap+WINDOWSIZE)));
    A_RMS(m,1) = (1/(WINDOWSIZE)) * sqrt(sum(y.*y)); % Normalization
    A_E(m,1) = (1/(WINDOWSIZE)) * sum(abs(y).^2); % Normalization
end

%% Calculating Zero-Crossing Rate for each window
% A_ZCR = zeros(nofw,1);
% for m = 1:nofw
%     zc = 0;
%     y = s1(round((m-1)*noverlap+1):min(length(s1), floor((m-1)*noverlap+WINDOWSIZE)));
%     for i=2:length(y)
%         if y(i)*y(i-1) < 0
%         zc = zc+1;
%         end
%     end
%     ZCR(m,1) = sum(abs(diff(y>0)))/length(y);
% end

%% Calculating Spectral Flux
H = hamming(WINDOWSIZE);
A_SF = zeros(nofw,1);
for m=1:nofw
    window = H.*(s1(round((m-1)*noverlap+1):min(length(s1), floor((m-1)*noverlap+WINDOWSIZE))));    
    FFT = (abs(fft(window,2*WINDOWSIZE)));
    FFT = FFT(1:WINDOWSIZE);        
    FFT = FFT / max(FFT);
    if (m>1)
        A_SF(m) = (1/(WINDOWSIZE)) * sum((FFT-FFTprev).^2); % Normalization
    else
        A_SF(m) = 0;
    end
    FFTprev = FFT;
end

%% Calculating Energy Entropy Block
% Eol = sum(s1.^2);
% numOfShortBlocks = 2;
% EnEntropy = zeros(nofw,1);
% curPos = 1;
% for m=1:nofw
%     curBlock = s1(curPos:curPos+WINDOWSIZE-1);
%     for j=1:numOfShortBlocks        
%         s(j) = sum(curBlock((j-1)*(WINDOWSIZE/numOfShortBlocks)+1:j*(WINDOWSIZE/numOfShortBlocks)).^2)/Eol;
%     end
%     
%     EnEntropy(m) = (1/(WINDOWSIZE)) * (-sum(s.*log2(s))); % Normalization
%     curPos = curPos + noverlap;
% end

%% Calculating Spectral Entropy

curPos = 1;
H = hamming(WINDOWSIZE);
A_SpecEntropy = zeros(nofw,1);
numOfBins = FFTSIZE;
h_step = FFTSIZE / numOfBins;
for i=1:nofw
    window = (H.*s1(curPos:curPos+WINDOWSIZE-1));
    fftTemp = abs(fft(window,2*FFTSIZE));
    fftTemp = fftTemp(1:FFTSIZE);
    S = sum(fftTemp);    
    
    for j=1:numOfBins
        x(j) = sum(fftTemp((j-1)*h_step + 1: j*h_step)) / S;
    end
    A_SpecEntropy(i) = (1/(WINDOWSIZE)) * (-sum(x.*log2(x)));
    curPos = curPos + noverlap;
end

%% Calculating Spectral Centroid

curPos = 1;
H = hamming(WINDOWSIZE);
m = ((fs/(2*WINDOWSIZE))*(1:WINDOWSIZE))';
A_SC = zeros(nofw,1);
for i=1:nofw
    window = H.*(s1(curPos:curPos+WINDOWSIZE-1));    
    FFT = (abs(fft(window,2*WINDOWSIZE)));
    FFT = FFT(1:WINDOWSIZE);  
    FFT = FFT / max(FFT);
    A_SC(i) = sum(m.*FFT)/sum(FFT);
    curPos = curPos + noverlap;
end
A_SC = A_SC / (fs/2); % Normalizing

%% Calculating Spectral Roll Off
c = 0.85;
curPos = 1;
H = hamming(WINDOWSIZE);
for i=1:nofw
    window = H.*(s1(curPos:curPos+WINDOWSIZE-1));    
    FFT = (abs(fft(window,512)));
    FFT = FFT(1:255);
    totalEnergy = sum(FFT);
    curEnergy = 0.0;
    countFFT = 1;
    while ((curEnergy<=c*totalEnergy) && (countFFT<=255))
        curEnergy = curEnergy + FFT(countFFT);
        countFFT = countFFT + 1;
    end
   A_SRO(i, 1) = ((countFFT-1))/(fs/2);
    curPos = curPos + noverlap;
end
%% Mel-Frequency Coefficients
% coeffs = mfcc(s1, fs,'WindowLength',WINDOWSIZE,'OverlapLength',noverlap, 'NumCoeffs',40, 'FFTLength',FFTSIZE);
% A_mean_coeffs = mean(coeffs, 2);
end
