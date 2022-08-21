% -----------------------------------------------------------------------
% Thesis
%
% Contact: Behnam Sherafat @ < behnam.sherafat@utah.edu>
% -----------------------------------------------------------------------
%% Initializing

clear;
close all;
clc;

% Markov Parameters

markovV_m = [0.990196078,	0.009803922; 0.017676768,	0.982323232];
sum(markovV_m,2)

markovA_m = [0.990196078,	0.009803922; 0.017676768,	0.982323232];
sum(markovA_m,2)

markovF_m = (markovV_m + markovA_m)/2;
sum(markovF_m,2)

% PCA dimension reduction factor
PCA_f = 0.5;

% Parameters for smoothing windows

szV_s = 2;  % size for small window
szV_b = 6;  % size for big window
threV = 0.6; % threshold

szA_s = 2;  % size for small window
szA_b = 6;  % size for big window
threA = 0.6; % threshold

szF_s = 2;  % size for small window
szF_b = 6;  % size for big window
threF = 0.6; % threshold

%% Reading Sensor Data

load test2.mat

fsAccel = length(timestampAccel)/timestampAccel(length(timestampAccel));
fsAng = length(timestampAng)/timestampAng(length(timestampAng));
fsMag = length(timestampMag)/timestampMag(length(timestampMag));
fsOri = length(timestampOri)/timestampOri(length(timestampOri));
fs = [fsAccel fsAng fsMag fsOri];
% fsV = floor(min(fs));
fsV = fsAccel;

tAccel = timestampAccel(length(timestampAccel));
tAng = timestampAng(length(timestampAng));
tMag = timestampMag(length(timestampMag));
tOri = timestampOri(length(timestampOri));
T1 = [tAccel tAng tMag tOri];
% t_V = min(T1); % in seconds

%% Removing Sensor Data outliers

% Removing G from Accelerometer Data
% Read Only Acceleration Sensor Data
X = Ax;
Y = Ay;
Z = Az;

t = timestampAccel;

kFilteringFactor = 0.1;
g1 = 0;
g2 = 0;
g3 = 0;
Bx = zeros(size(X, 1), 1);
By = zeros(size(Y, 1), 1);
Bz = zeros(size(Z, 1), 1);

for i = 1:size(X, 1)
    
    g1 = (X(i, 1) * kFilteringFactor) + (g1 * (1.0 - kFilteringFactor));
    g2 = (Y(i, 1) * kFilteringFactor) + (g2 * (1.0 - kFilteringFactor));
    g3 = (Z(i, 1) * kFilteringFactor) + (g3 * (1.0 - kFilteringFactor));

    Bx(i, 1) = X(i, 1) - g1;
    By(i, 1) = Y(i, 1) - g2;
    Bz(i, 1) = Z(i, 1) - g3;
end

[Refined_Ax,~,~,~,~] = filloutliers(Bx,'clip','movmedian', 3,'SamplePoints',t);
[Refined_Ay,~,~,~,~] = filloutliers(By,'clip','movmedian', 3,'SamplePoints',t);
[Refined_Az,~,~,~,~] = filloutliers(Bz,'clip','movmedian', 3,'SamplePoints',t);

Refilled_Ax = fillmissing(Refined_Ax,'movmedian',24); 
Refilled_Ay = fillmissing(Refined_Ay,'movmedian',24); 
Refilled_Az = fillmissing(Refined_Az,'movmedian',24); 

MagAccel = sqrt(sum(Refilled_Ax.^2 + Refilled_Ay.^2 + Refilled_Az.^2, 2));
MagAccel = MagAccel - mean(MagAccel);
Mag_A = MagAccel/max(MagAccel);

S_Vibration = Refined_Ax; % Other options: Refined_Ax or Refined_Ay or Refined_Az or Mag_A

% Read other Sensor Data
% signal=MagAng;
% t = timestampAng;
% fs = length(signal)/t(length(t));
% signal = signal - mean(signal);
% Mag_A = signal/max(signal);
% final_Signal = Mag_A;

%% Reading Audio Data

[A_ori, fsA] = audioread('ZOOM0030.wav');

%% Audio Denoising: method = 'mcra2'; 

A_ori(A_ori==0)=10^-4;
A_ori = A_ori(:, 1);
audiowrite('ZOOM0030_ori.wav', A_ori, fsA);

% Audio Denoising: method = 'mcra2'; 
% specsub_ns('ZOOM0030_ori.wav', 'mcra2', 'ZOOM0030_den.wav')

% Read Modified Audio

[S_Audio, fsA] = audioread('ZOOM0030_den.wav');
% t_A = size(S_Audio, 1)/fsA;


%% Mutual Time Stamp and Sycnhronization

% start_Audio = round(0.38*fsA);
% S_Audio = S_Audio(start_Audio:end);
% t_A = size(S_Audio, 1)/fsA;
% S_Audio = S_Audio - mean(S_Audio);
% start_Audio = 0;
% 
start_Vibration = round(0.1*fsV);
S_Vibration = S_Vibration(start_Vibration:end);
t_V = size(S_Vibration, 1)/fsV;
start_Vibration = 0;

% 
% start_Vibration = 0;
% S_Vibration = S_Vibration(start_Vibration+1:end);
% t_V = size(S_Vibration, 1)/fsV;
% 
start_Audio = 0;
S_Audio = S_Audio(start_Audio+1:end);
t_A = size(S_Audio, 1)/fsA;

T2 = [t_V t_A];
t_F = min(floor(T2)); % fusion mutual time stamp

S_Audio = S_Audio(start_Audio+1:floor(t_F*fsA));
S_Vibration = S_Vibration(start_Vibration+1:floor(t_F*fsV));


%% Defining STFT Parameters

H = size(S_Audio, 1)/size(S_Vibration, 1);

V_noverlap = 8;
V_WINDOWSIZE = V_noverlap * 2;

A_noverlap = round(size(S_Audio, 1)/(1+((size(S_Vibration, 1)-V_noverlap)/V_noverlap)));
% A_noverlap = ceil(4*H);
A_WINDOWSIZE = A_noverlap * 2 ;

FFTSIZE = 48;

%% Plotting Data

figure(1);
subplot(4,1,1)
plot(S_Vibration)
xt1 = get(gca, 'XTick');                                 % 'XTick' Values
set(gca,'YDir','normal');
set(gca, 'XTick', xt1, 'XTickLabel', round(xt1/fsV))
title('CAT 259D Vibration Signal','FontSize',20,'FontWeight','bold');
y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight


subplot(4,1,2)
[V_Spec, ~, ~, ~, ~, ~, ~, ~] = FunctionExtractionV(S_Vibration, 0, t_V, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
imagesc(db(V_Spec')); colormap('jet'); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 259D Vibration Signal Spectrogram',...
    'FontSize', 20, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
axis tight


subplot(4,1,3)
plot(S_Audio(1:74*fsA))
xt2 = get(gca, 'XTick');                                 % 'XTick' Values
set(gca,'YDir','normal');
set(gca, 'XTick', xt2, 'XTickLabel', round(xt2/fsA))
title('CAT 259D Audio Signal','FontSize',20,'FontWeight','bold');
y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight

subplot(4,1,4)
[A_Spec, ~, ~, ~, ~, ~, ~, ~] = FunctionExtractionA(S_Audio, 0, t_A, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
imagesc(db(A_Spec')); colormap('jet'); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 259D Audio Signal Spectrogram',...
    'FontSize', 20, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
axis tight

%% Building Vibration and Audio Data Set
% ================================ Building Major Library  ===========================
s_start = 7; % start time 
s_end = 8; % end time
% ============================== Training Vibration 1 ============================== %

[V_Spec, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vibration, s_start, s_end, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
F_V1 = [V_Spec, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

figure(2); 
subplot(5, 1, 1);
imagesc(db(V_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram for Vibration Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s'], 'FontSize',20,'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight

% ============================== Training Audio 1 ============================== %

[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, s_start, s_end, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
F_A1 = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];

figure(4); 
subplot(5, 1, 1);
imagesc(db(A_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram for Audio Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s'], 'FontSize',20,'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight



s_start = 34; % start time
s_end = 41; % end time

% ============================== Training Vibration 2 ============================== %
[V_Spec, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vibration, s_start, s_end, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
F_V2 = [V_Spec, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

figure(2); 
subplot(5, 1, 2);
imagesc(db(V_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');

% ============================== Training Audio 2 ============================== %

[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, s_start, s_end, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
F_A2 = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];

figure(4); 
subplot(5, 1, 2);
imagesc(db(A_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');


s_start = 50.6; % start time
s_end = 51.4; % end time

% ============================== Training Vibration 3 ============================== %

[V_Spec, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vibration, s_start, s_end, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
F_V3 = [V_Spec, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

figure(2); 
subplot(5, 1, 3);
imagesc(db(V_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');

% ============================== Training Audio 3 ============================== %

[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, s_start, s_end, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
F_A3 = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];

figure(4); 
subplot(5, 1, 3);
imagesc(db(A_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');



s_start = 60.5; % start time
s_end = 62; % end time

% ============================== Training Vibration 4 ============================== %

[V_Spec, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vibration, s_start, s_end, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
F_V4 = [V_Spec, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];


figure(2); 
subplot(5, 1, 4);
imagesc(db(V_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');

% ============================== Training Audio 4 ============================== %

[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, s_start, s_end, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
F_A4 = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];


figure(4); 
subplot(5, 1, 4);
imagesc(db(A_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');

% ==================================================================== %

S_concaV1 = [abs(F_V1); abs(F_V2); abs(F_V3); abs(F_V4)]; % concatenate all STFT matrix

figure(2); 
subplot(5, 1, 5)
imagesc(db(S_concaV1')); set(gca,'YDir','normal', 'Clim', c);
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Spectrogram Major: Concatenate above')
ylabel('Normalized Frequency'); xlabel('Time');

% ==================================================================== %

S_concaA1 = [abs(F_A1); abs(F_A2); abs(F_A3); abs(F_A4)]; % concatenate all STFT matrix

figure(4); 
subplot(5, 1, 5)
imagesc(db(S_concaA1')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Spectrogram Major: Concatenate above')
ylabel('Normalized Frequency'); xlabel('Time');
 
% =================  Building Minor Library ===========================

s_start = 0; % start time 
s_end = 5; % end time

% ============================== Training Vibration 5 ============================== %

[V_Spec, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vibration, s_start, s_end, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
F_V5 = [V_Spec, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

figure(3); 
subplot(5, 1, 1);
imagesc(db(V_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');

% ============================== Training Audio 5 ============================== %

[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, s_start, s_end, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
F_A5 = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];


figure(5); 
subplot(5, 1, 1);
imagesc(db(A_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');


s_start = 46; % start time
s_end = 50; % end time

% ============================== Training Vibration 6 ============================== %


[V_Spec, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vibration, s_start, s_end, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
F_V6 = [V_Spec, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];
figure(3); 
subplot(5, 1, 2);
imagesc(db(V_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');

% ============================== Training Audio 6 ============================== %

[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, s_start, s_end, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
F_A6 = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];


figure(5); 
subplot(5, 1, 2);
imagesc(db(A_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');


s_start = 52; % start time
s_end = 53; % end time

% ============================== Training Vibration 7 ============================== %

[V_Spec, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vibration, s_start, s_end, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
F_V7 = [V_Spec, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

figure(3); 
subplot(5, 1, 3);
imagesc(db(V_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');

% ============================== Training Audio 7 ============================== %

[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, s_start, s_end, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
F_A7 = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];

figure(5); 
subplot(5, 1, 3);
imagesc(db(A_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');


s_start = 55; % start time
s_end = 57; % end time

% ============================== Training Vibration 8 ============================== %


[V_Spec, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vibration, s_start, s_end, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
F_V8 = [V_Spec, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

figure(3); 
subplot(5, 1, 4);
imagesc(db(V_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');

% ============================== Training Audio 8 ============================== %

[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, s_start, s_end, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
F_A8 = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];

figure(5); 
subplot(5, 1, 4);
imagesc(db(A_Spec')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title(['Spectrogram Major Sample 1: ', num2str(s_start), '-', num2str(s_end), ' s']);
ylabel('Normalized Frequency'); xlabel('Time');

% ==================================================================== %

S_concaV2 = [abs(F_V5); abs(F_V6); abs(F_V7); abs(F_V8)]; % concatenate all STFT matrix

figure(3); 
subplot(5, 1, 5)
imagesc(db(S_concaV2')); set(gca,'YDir','normal', 'Clim', c);
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Spectrogram Minor: Concatenate above')
ylabel('Normalized Frequency'); xlabel('Time');

% ==================================================================== %

S_concaA2 = [abs(F_A5); abs(F_A6); abs(F_A7); abs(F_A8)]; % concatenate all STFT matrix

figure(5); 
subplot(5, 1, 5)
imagesc(db(S_concaA2')); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('Spectrogram Minor: Concatenate above')
ylabel('Normalized Frequency'); xlabel('Time');

%% Check Features Seperability
% 
% Data1 = S_concaA1;
% Data2 = S_concaA2;
% [H1, H2] = computeHistError(Data1, Data2);

%% Check Vibration and Audio Observations Equality

S_TrainV = [S_concaV1; S_concaV2];
S_TrainA = [S_concaA1; S_concaA2];

if size(S_TrainV,1) == size(S_TrainA,1)
    N_Observations = size(S_TrainV,1);
else
    warning('Number of observations in audio and vibratio do not match');
end

%% Building Vibration Data Test

% TrainingV_Mean = mean([S_concaV1; S_concaV2], 1);
% S_concaV1 = (S_concaV1-TrainingV_Mean);
% S_concaV2 = (S_concaV2-TrainingV_Mean);
% 
S_majorV = S_concaV1;
S_minorV = S_concaV2;

l_LibV1 = 2.*ones(size(S_concaV1, 1), 1);
l_LibV2 = ones(size(S_concaV2, 1), 1);

TrainingSampleV = [S_majorV; S_minorV];
TrainingLabelV = [l_LibV1; l_LibV2];

% Implementing PCA
p = round(0.8*size(TrainingSampleV, 2));
[COEFF, SCORE, pcvar] = pca(TrainingSampleV);
p = round(PCA_f*size(TrainingSampleV, 2));
Ap = COEFF(:, 1:p);
RTrainingSampleV = bsxfun(@minus, TrainingSampleV, mean(TrainingSampleV))*Ap;

% SVM Training
SVMModelV = fitcsvm(RTrainingSampleV,TrainingLabelV,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto');

% start and finish time of test file
s_test_startV = 0;

% Extract Features
[V_Spec, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vibration, start_Vibration, t_F, fsV, FFTSIZE, V_noverlap, V_WINDOWSIZE);
F_TV = [V_Spec, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

% Add Actual Labels

label_act = 2.*ones(size(V_Spec', 2), 1);
lps = length(label_act)/(74-start_Vibration); % length per second 

minor_start1 = ceil(start_Vibration*lps+1);
minor_end1 = ceil(6.1*lps);
minor_start2 = ceil(8.5*lps);
minor_end2 = ceil(10.5*lps);
minor_start3 = ceil(13.4*lps);
minor_end3 = ceil(15.45*lps);
minor_start4 = ceil(17.1*lps);
minor_end4 = ceil(19.4*lps);
minor_start5 = ceil(45*lps);
minor_end5 = ceil(50.55*lps);
minor_start6 = ceil(51.5*lps);
minor_end6 = ceil(53.85*lps);
minor_start7 = ceil(54.65*lps);
minor_end7 = ceil(57*lps);
minor_start8 = ceil(58.25*lps);
minor_end8 = ceil(60.2*lps);
minor_start9 = ceil(62*lps);
minor_end9 = ceil(66.2*lps-1);

label_act(minor_start1:minor_end1) = 1;
label_act(minor_start2:minor_end2) = 1;
label_act(minor_start3:minor_end3) = 1;
label_act(minor_start4:minor_end4) = 1;
label_act(minor_start5:minor_end5) = 1;
label_act(minor_start6:minor_end6) = 1;
label_act(minor_start7:minor_end7) = 1;
label_act(minor_start8:minor_end8) = 1;
label_act(minor_start9:minor_end9) = 1;


% SVM Testing
% F_TV = (F_TV-TrainingV_Mean);
TestingSampleV = F_TV;
RTestingSampleV = bsxfun(@minus, TestingSampleV, mean(TrainingSampleV))*Ap;
[predictedV_label] = predict(SVMModelV, RTestingSampleV);
x = 1:length(predictedV_label);

% ================================= Small window ================================= %
filteredV_sW = predictedV_label;

labelsV = unique(predictedV_label);   % Find how many different labels
nlabelsV = zeros(1, length(labelsV)); % matrix used to record numbers of each label
accV_sw=zeros(length(predictedV_label),1);


for ii = 1:length(predictedV_label)
    % Beginning
    if ii <= szV_s
        for jj = 1:length(labelsV)
            nlabelsV(jj) = length(find(predictedV_label(1:ii+szV_s) == labelsV(jj)));
            if(nlabelsV(jj)/(ii+szV_s) >= threV)
                filteredV_sW(ii) = labelsV(jj);
                accV_sw(ii)=nlabelsV(jj)/(ii+szV_s); 
            end
        end
    % Ending
    elseif ii >= length(filteredV_sW) - szV_s
        for jj = 1:length(labelsV)
            nlabelsV(jj) = length(find(predictedV_label(ii-szV_s:end) == labelsV(jj)));
            if(nlabelsV(jj)/(length(predictedV_label)+szV_s-ii+1) >= threV)
                filteredV_sW(ii) = labelsV(jj);
                 accV_sw(ii)=nlabelsV(jj)/(length(predictedV_label)+szV_s-ii+1);
            end
        end
    % Mid
    else
        for jj = 1:length(labelsV)
            nlabelsV(jj) = length(find(predictedV_label(ii-szV_s:ii+szV_s) == labelsV(jj)));
            if(nlabelsV(jj)/(2*szV_s+1) >= threV)
                filteredV_sW(ii) = labelsV(jj);
                accV_sw(ii)=nlabelsV(jj)/(2*szV_s+1); 
            end
        end
    end
end

% ================================= Large window ================================= %
filteredV_bW = filteredV_sW; % Initilize labels filtered by large window
accV_bw=zeros(length(filteredV_bW),1);
% Large window
for ii = 1:length(filteredV_bW)
    % Beginning
    if ii <= szV_b
        for jj = 1:length(labelsV)
            nlabelsV(jj) = length(find(filteredV_sW(1:ii+szV_b) == labelsV(jj)));
            if(nlabelsV(jj)/(ii+szV_b) >= threV)
                filteredV_bW(ii) = labelsV(jj);
                accV_bw(ii)=(nlabelsV(jj)/(ii+szV_b));
            end
        end
    % Ending
    elseif ii >= length(filteredV_bW) - szV_b
        for jj = 1:length(labelsV)
            nlabelsV(jj) = length(find(filteredV_sW(ii-szV_b:end) == labelsV(jj)));
            if(nlabelsV(jj)/(length(filteredV_sW)+szV_b-ii+1) >= threV)
                filteredV_bW(ii) = labelsV(jj);
                accV_bw(ii)=(nlabelsV(jj)/(length(filteredV_sW)+szV_b-ii+1));
            end
        end
    % Mid
    else
        for jj = 1:length(labelsV)
            nlabelsV(jj) = length(find(filteredV_sW(ii-szV_b:ii+szV_b) == labelsV(jj)));
            if(nlabelsV(jj)/(2*szV_b+1) >= threV)
                filteredV_bW(ii) = labelsV(jj);
                accV_bw(ii)=(nlabelsV(jj)/(2*szV_b+1));
            end
        end
    end
end

% Performance matrix Audio
pV1 = zeros(length(labelsV));

% ===================================== %
% Performance matrix:  p1
%
%                   correct label
%                     1        2
%                 -----------------
%              1  - 1/1  -  2/1   -sum
%   predicted     ----------------- 
%              2  - 1/2  -  2/2   -
%                 -----------------
% ===================================== %

for ii = 1:length(label_act)
    if(filteredV_bW(ii) == 1 && label_act(ii) == 1)
        pV1(1,1) = pV1(1,1) +1;
    elseif(filteredV_bW(ii) == 1 && label_act(ii) == 2)
        pV1(1,2) = pV1(1,2) +1;
    elseif(filteredV_bW(ii) == 2 && label_act(ii) == 1)
        pV1(2,1) = pV1(2,1) +1;
    else
        pV1(2,2) = pV1(2,2) +1;
    end
end

pV_1(1,1) = pV1(1,1)/(pV1(1,1)+pV1(2,1));
pV_1(2,1) = pV1(2,1)/(pV1(1,1)+pV1(2,1));
pV_1(1,2) = pV1(1,2)/(pV1(1,2)+pV1(2,2));
pV_1(2,2) = pV1(2,2)/(pV1(1,2)+pV1(2,2));

save('PerformanceMatrixV1.mat', 'pV1', '-mat');

% Markov

% preallocation
outstateV=filteredV_sW;
accV_mk=zeros(length(accV_sw),1);

% start period
T=1;
C=0;

for ii= 3:length(filteredV_bW)
    
    prev_state=(outstateV(ii-1));
    pred_state=(filteredV_bW(ii));
    
  if prev_state==1
        prev_m=[1,0];
    else
        prev_m=[0,1];
  end
    
  if outstateV(ii-1)==outstateV(ii-2)
        T=T+1;
    else
        T=1;
        C=C+1;
  end

    markovV=prev_m*(markovV_m^T);
if markovV(1,1)>0.5 
    accV_mk(ii)=markovV(1,1);
    markovV=[1,0];
    p_markov=1;
else
    if markovV(1,2)>0.5
        accV_mk(ii)=markovV(1,2);
        markovV=[0,1];
        p_markov=2;
    else
        p_markov=pred_state;
    end
end

%Decision.

if accV_mk(ii)>=accV_bw(ii)
    outstateV(ii)=p_markov;
else
    outstateV(ii)=pred_state;
end
end

%
figure(7); 
s(1) = subplot(4, 1, 1);
imagesc(db(V_Spec')); colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Vibration "Predicted SW"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*filteredV_sW, 'filled', 'Marker', 'none');
% plot(x, 100.*filtered_bW, 'ob', x, 100.*label_ans, 'sk');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');

s(2) = subplot(4, 1, 2);
imagesc(db(V_Spec')); colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Vibration "Predicted BW"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*filteredV_bW, 'filled', 'Marker', 'none');
% plot(x, 100.*filtered_bW, 'ob', x, 100.*label_ans, 'sk');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');


s(3) = subplot(4, 1, 3);
imagesc(db(V_Spec')); colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Vibration "Predicted Markov"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*outstateV, 'filled', 'Marker', 'none');
% plot(x, 100.*filtered_bW, 'ob', x, 100.*label_ans, 'sk');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');

s(4) = subplot(4, 1, 4); 
imagesc(db(V_Spec'));% colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Vibration "Correct Label"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
% stem(x, 100.*label_ans, 'filled', 'Marker', 'none');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');

%
pV2 = zeros(length(labelsV));

% ===================================== %
% Performance matrix:  p
%
%                   correct label
%                     1        2
%                 -----------------
%              1  - 1/1  -  2/1   -sum
%   predicted     ----------------- 
%              2  - 1/2  -  2/2   -
%                 -----------------
% ===================================== %

for ii = 1:length(label_act)
    if(outstateV(ii) == 1 && label_act(ii) == 1)
        pV2(1,1) = pV2(1,1) +1;
    elseif(outstateV(ii) == 1 && label_act(ii) == 2)
        pV2(1,2) = pV2(1,2) +1;
    elseif(outstateV(ii) == 2 && label_act(ii) == 1)
        pV2(2,1) = pV2(2,1) +1;
    else
        pV2(2,2) = pV2(2,2) +1;
    end
end

pV_2(1,1) = pV2(1,1)/(pV2(1,1)+pV2(2,1));
pV_2(2,1) = pV2(2,1)/(pV2(1,1)+pV2(2,1));
pV_2(1,2) = pV2(1,2)/(pV2(1,2)+pV2(2,2));
pV_2(2,2) = pV2(2,2)/(pV2(1,2)+pV2(2,2));

save('PerformanceMatrixV2.mat', 'pV2', '-mat');


fprintf('After Window Filtering \n');
pV1
fprintf('After Markov Chain \n');
pV2
AccuracyV = (pV2(1,1)+pV2(2,2))/((pV2(1,1)+pV2(1,2)+(pV2(2,1)+pV2(2,2))))
%% Building Audio Data Test

% TrainingA_Mean = mean([S_concaA1; S_concaA2], 1);
% S_concaA1 = (S_concaA1-TrainingA_Mean);
% S_concaA2 = (S_concaA2-TrainingA_Mean);

S_majorA = S_concaA1;
S_minorA = S_concaA2;

l_LibA1 = 2.*ones(size(S_concaA1, 1), 1);
l_LibA2 = ones(size(S_concaA2, 1), 1);

TrainingSampleA = [S_majorA; S_minorA];
TrainingLabelA = [l_LibA1; l_LibA2];

% Implementing PCA

[COEFF, SCORE] = pca(TrainingSampleA);
p = round(PCA_f*size(TrainingSampleA, 2));
Ap = COEFF(:, 1:p);
RTrainingSampleA = bsxfun(@minus, TrainingSampleA, mean(TrainingSampleA))*Ap;

% SVM Training
SVMModelA = fitcsvm(RTrainingSampleA,TrainingLabelA,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto');

% Extract Features
[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, start_Audio, t_F, fsA, FFTSIZE, A_noverlap, A_WINDOWSIZE);
F_TA = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];

% SVM Testing
% F_TA = (F_TA-TrainingA_Mean);
TestingSampleA = F_TA;
RTestingSampleA = bsxfun(@minus, TestingSampleA, mean(TrainingSampleA))*Ap;
[predictedA_label] = predict(SVMModelA, RTestingSampleA);
x = 1:length(predictedA_label);

% ================================= Small window ================================= %
filteredA_sW = predictedA_label;

labelsA = unique(predictedA_label);   % Find how many different labels
nlabelsA = zeros(1, length(labelsA)); % matrix used to record numbers of each label
accA_sw=zeros(length(predictedA_label),1);


for ii = 1:length(predictedA_label)
    % Beginning
    if ii <= szA_s
        for jj = 1:length(labelsA)
            nlabelsA(jj) = length(find(predictedA_label(1:ii+szA_s) == labelsA(jj)));
            if(nlabelsA(jj)/(ii+szA_s) >= threA)
                filteredA_sW(ii) = labelsA(jj);
                accA_sw(ii)=nlabelsA(jj)/(ii+szA_s); 
            end
        end
    % Ending
    elseif ii >= length(filteredA_sW) - szA_s
        for jj = 1:length(labelsA)
            nlabelsA(jj) = length(find(predictedA_label(ii-szA_s:end) == labelsA(jj)));
            if(nlabelsA(jj)/(length(predictedA_label)+szA_s-ii+1) >= threA)
                filteredA_sW(ii) = labelsA(jj);
                 accA_sw(ii)=nlabelsA(jj)/(length(predictedA_label)+szA_s-ii+1);
            end
        end
    % Mid
    else
        for jj = 1:length(labelsA)
            nlabelsA(jj) = length(find(predictedA_label(ii-szA_s:ii+szA_s) == labelsA(jj)));
            if(nlabelsA(jj)/(2*szA_s+1) >= threA)
                filteredA_sW(ii) = labelsA(jj);
                accA_sw(ii)=nlabelsA(jj)/(2*szA_s+1); 
            end
        end
    end
end

% ================================= Large window ================================= %
filteredA_bW = filteredA_sW; % Initilize labels filtered by large window
accA_bw=zeros(length(filteredA_bW),1);
% Large window
for ii = 1:length(filteredA_bW)
    % Beginning
    if ii <= szA_b
        for jj = 1:length(labelsA)
            nlabelsA(jj) = length(find(filteredA_sW(1:ii+szA_b) == labelsA(jj)));
            if(nlabelsA(jj)/(ii+szA_b) >= threA)
                filteredA_bW(ii) = labelsA(jj);
                accA_bw(ii)=(nlabelsA(jj)/(ii+szA_b));
            end
        end
    % Ending
    elseif ii >= length(filteredA_bW) - szA_b
        for jj = 1:length(labelsA)
            nlabelsA(jj) = length(find(filteredA_sW(ii-szA_b:end) == labelsA(jj)));
            if(nlabelsA(jj)/(length(filteredA_sW)+szA_b-ii+1) >= threA)
                filteredA_bW(ii) = labelsA(jj);
                accA_bw(ii)=(nlabelsA(jj)/(length(filteredA_sW)+szA_b-ii+1));
            end
        end
    % Mid
    else
        for jj = 1:length(labelsA)
            nlabelsA(jj) = length(find(filteredA_sW(ii-szA_b:ii+szA_b) == labelsA(jj)));
            if(nlabelsA(jj)/(2*szA_b+1) >= threA)
                filteredA_bW(ii) = labelsA(jj);
                accA_bw(ii)=(nlabelsA(jj)/(2*szA_b+1));
            end
        end
    end
end

% Performance matrix Audio
pA1 = zeros(length(labelsA));

% ===================================== %
% Performance matrix:  p1
%
%                   correct label
%                     1        2
%                 -----------------
%              1  - 1/1  -  2/1   -sum
%   predicted     ----------------- 
%              2  - 1/2  -  2/2   -
%                 -----------------
% ===================================== %

for ii = 1:length(label_act)
    if(filteredA_bW(ii) == 1 && label_act(ii) == 1)
        pA1(1,1) = pA1(1,1) +1;
    elseif(filteredA_bW(ii) == 1 && label_act(ii) == 2)
        pA1(1,2) = pA1(1,2) +1;
    elseif(filteredA_bW(ii) == 2 && label_act(ii) == 1)
        pA1(2,1) = pA1(2,1) +1;
    else
        pA1(2,2) = pA1(2,2) +1;
    end
end

pA_1(1,1) = pA1(1,1)/(pA1(1,1)+pA1(2,1));
pA_1(2,1) = pA1(2,1)/(pA1(1,1)+pA1(2,1));
pA_1(1,2) = pA1(1,2)/(pA1(1,2)+pA1(2,2));
pA_1(2,2) = pA1(2,2)/(pA1(1,2)+pA1(2,2));

save('PerformanceMatrixA1.mat', 'pA1', '-mat');

% Markov

% preallocation
outstateA=filteredA_sW;
accA_mk=zeros(length(accA_sw),1);

% start period
T=1;
C=0;

for ii= 3:length(filteredA_bW)
    
    prev_state=(outstateA(ii-1));
    pred_state=(filteredA_bW(ii));
    
  if prev_state==1
        prev_m=[1,0];
    else
        prev_m=[0,1];
  end
    
  if outstateA(ii-1)==outstateA(ii-2)
        T=T+1;
    else
        T=1;
        C=C+1;
  end

    markovA=prev_m*(markovA_m^T);
if markovA(1,1)>0.5 
    accA_mk(ii)=markovA(1,1);
    markovA=[1,0];
    p_markov=1;
else
    if markovA(1,2)>0.5
        accA_mk(ii)=markovA(1,2);
        markovA=[0,1];
        p_markov=2;
    else
        p_markov=pred_state;
    end
end

%Decision.

if accA_mk(ii)>=accA_bw(ii)
    outstateA(ii)=p_markov;
else
    outstateA(ii)=pred_state;
end
end

%%
figure(6); 
s(1) = subplot(4, 1, 1);
imagesc(db(A_Spec')); colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Audio "Predicted SW"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*filteredA_sW, 'filled', 'Marker', 'none');
% plot(x, 100.*filtered_bW, 'ob', x, 100.*label_ans, 'sk');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');

s(2) = subplot(4, 1, 2);
imagesc(db(A_Spec')); colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Audio "Predicted BW"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*filteredA_bW, 'filled', 'Marker', 'none');
% plot(x, 100.*filtered_bW, 'ob', x, 100.*label_ans, 'sk');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');


s(3) = subplot(4, 1, 3);
imagesc(db(A_Spec')); colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Audio "Predicted Markov"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*outstateA, 'filled', 'Marker', 'none');
% plot(x, 100.*filtered_bW, 'ob', x, 100.*label_ans, 'sk');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');

s(4) = subplot(4, 1, 4); 
imagesc(db(A_Spec'));% colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Audio "Correct Label"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
% stem(x, 100.*label_ans, 'filled', 'Marker', 'none');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');

%
pA2 = zeros(length(labelsA));

% ===================================== %
% Performance matrix:  p
%
%                   correct label
%                     1        2
%                 -----------------
%              1  - 1/1  -  2/1   -sum
%   predicted     ----------------- 
%              2  - 1/2  -  2/2   -
%                 -----------------
% ===================================== %

for ii = 1:length(label_act)
    if(outstateA(ii) == 1 && label_act(ii) == 1)
        pA2(1,1) = pA2(1,1) +1;
    elseif(outstateA(ii) == 1 && label_act(ii) == 2)
        pA2(1,2) = pA2(1,2) +1;
    elseif(outstateA(ii) == 2 && label_act(ii) == 1)
        pA2(2,1) = pA2(2,1) +1;
    else
        pA2(2,2) = pA2(2,2) +1;
    end
end

pA_2(1,1) = pA2(1,1)/(pA2(1,1)+pA2(2,1));
pA_2(2,1) = pA2(2,1)/(pA2(1,1)+pA2(2,1));
pA_2(1,2) = pA2(1,2)/(pA2(1,2)+pA2(2,2));
pA_2(2,2) = pA2(2,2)/(pA2(1,2)+pA2(2,2));

save('PerformanceMatrixA2.mat', 'pA2', '-mat');


fprintf('After Window Filtering \n');
pA1
fprintf('After Markov Chain \n');
pA2
AccuracyA = (pA2(1,1)+pA2(2,2))/((pA2(1,1)+pA2(1,2)+(pA2(2,1)+pA2(2,2))))
%% Building Fusion Data Test

% TrainingF_Mean = mean([S_concaA1, S_concaV1; S_concaA2, S_concaV2], 1);

S_concaF1 = [S_concaA1, S_concaV1];
S_concaF2 = [S_concaA2, S_concaV2];

% S_concaF1 = (S_concaF1-TrainingF_Mean);
% S_concaF2 = (S_concaF2-TrainingF_Mean);

S_majorF = S_concaF1;
S_minorF = S_concaF2;

l_LibF1 = 2.*ones(size(S_concaF1, 1), 1);
l_LibF2 = ones(size(S_concaF2, 1), 1);

TrainingSampleF = [S_majorF; S_minorF];
TrainingLabelF = [l_LibF1; l_LibF2];

% Implementing PCA

[COEFF, SCORE] = pca(TrainingSampleF);
p = round(PCA_f*size(TrainingSampleF, 2));
Ap = COEFF(:, 1:p);
RTrainingSampleF = bsxfun(@minus, TrainingSampleF, mean(TrainingSampleF))*Ap;

% SVM Training
SVMModelF = fitcsvm(RTrainingSampleF,TrainingLabelF,'Standardize',true,'KernelFunction','RBF',...
    'KernelScale','auto');

% start and finish time of test file
s_test_startF = 0;
s_test_endA = t_F;

% SVM Testing

F_TF = [F_TA, F_TV];
% F_TF = (F_TF-TrainingF_Mean);
TestingSampleF = F_TF;
RTestingSampleF = bsxfun(@minus, TestingSampleF, mean(TrainingSampleF))*Ap;

[predictedF_label] = predict(SVMModelF, RTestingSampleF);
x = 1:length(predictedF_label);

% ================================= Small window ================================= %
filteredF_sW = predictedF_label;

labelsF = unique(predictedF_label);   % Find how many different labels
nlabelsF = zeros(1, length(labelsF)); % matrix used to record numbers of each label
accF_sw=zeros(length(predictedF_label),1);


for ii = 1:length(predictedF_label)
    % Beginning
    if ii <= szF_s
        for jj = 1:length(labelsF)
            nlabelsF(jj) = length(find(predictedF_label(1:ii+szF_s) == labelsF(jj)));
            if(nlabelsF(jj)/(ii+szF_s) >= threF)
                filteredF_sW(ii) = labelsF(jj);
                accF_sw(ii)=nlabelsF(jj)/(ii+szF_s); 
            end
        end
    % Ending
    elseif ii >= length(filteredF_sW) - szF_s
        for jj = 1:length(labelsF)
            nlabelsF(jj) = length(find(predictedF_label(ii-szF_s:end) == labelsF(jj)));
            if(nlabelsF(jj)/(length(predictedF_label)+szF_s-ii+1) >= threF)
                filteredF_sW(ii) = labelsF(jj);
                 accF_sw(ii)=nlabelsF(jj)/(length(predictedF_label)+szF_s-ii+1);
            end
        end
    % Mid
    else
        for jj = 1:length(labelsF)
            nlabelsF(jj) = length(find(predictedF_label(ii-szF_s:ii+szF_s) == labelsF(jj)));
            if(nlabelsF(jj)/(2*szF_s+1) >= threF)
                filteredF_sW(ii) = labelsF(jj);
                accF_sw(ii)=nlabelsF(jj)/(2*szF_s+1); 
            end
        end
    end
end

% ================================= Large window ================================= %
filteredF_bW = filteredF_sW; % Initilize labels filtered by large window
accF_bw=zeros(length(filteredF_bW),1);
% Large window
for ii = 1:length(filteredF_bW)
    % Beginning
    if ii <= szF_b
        for jj = 1:length(labelsF)
            nlabelsF(jj) = length(find(filteredF_sW(1:ii+szF_b) == labelsF(jj)));
            if(nlabelsF(jj)/(ii+szF_b) >= threF)
                filteredF_bW(ii) = labelsF(jj);
                accF_bw(ii)=(nlabelsF(jj)/(ii+szF_b));
            end
        end
    % Ending
    elseif ii >= length(filteredF_bW) - szF_b
        for jj = 1:length(labelsF)
            nlabelsF(jj) = length(find(filteredF_sW(ii-szF_b:end) == labelsF(jj)));
            if(nlabelsF(jj)/(length(filteredF_sW)+szF_b-ii+1) >= threF)
                filteredF_bW(ii) = labelsF(jj);
                accF_bw(ii)=(nlabelsF(jj)/(length(filteredF_sW)+szF_b-ii+1));
            end
        end
    % Mid
    else
        for jj = 1:length(labelsF)
            nlabelsF(jj) = length(find(filteredF_sW(ii-szF_b:ii+szF_b) == labelsF(jj)));
            if(nlabelsF(jj)/(2*szF_b+1) >= threF)
                filteredF_bW(ii) = labelsF(jj);
                accF_bw(ii)=(nlabelsF(jj)/(2*szF_b+1));
            end
        end
    end
end

% Performance matrix Audio
pF1 = zeros(length(labelsF));

% ===================================== %
% Performance matrix:  p1
%
%                   correct label
%                     1        2
%                 -----------------
%              1  - 1/1  -  2/1   -sum
%   predicted     ----------------- 
%              2  - 1/2  -  2/2   -
%                 -----------------
% ===================================== %

for ii = 1:length(label_act)
    if(filteredF_bW(ii) == 1 && label_act(ii) == 1)
        pF1(1,1) = pF1(1,1) +1;
    elseif(filteredF_bW(ii) == 1 && label_act(ii) == 2)
        pF1(1,2) = pF1(1,2) +1;
    elseif(filteredF_bW(ii) == 2 && label_act(ii) == 1)
        pF1(2,1) = pF1(2,1) +1;
    else
        pF1(2,2) = pF1(2,2) +1;
    end
end

pF_1(1,1) = pF1(1,1)/(pF1(1,1)+pF1(2,1));
pF_1(2,1) = pF1(2,1)/(pF1(1,1)+pF1(2,1));
pF_1(1,2) = pF1(1,2)/(pF1(1,2)+pF1(2,2));
pF_1(2,2) = pF1(2,2)/(pF1(1,2)+pF1(2,2));

save('PerformanceMatrixF1.mat', 'pF1', '-mat');

% Markov

% preallocation
outstateF=filteredF_sW;
accF_mk=zeros(length(accF_sw),1);

% start period
T=1;
C=0;

for ii= 3:length(filteredF_bW)
    
    prev_state=(outstateF(ii-1));
    pred_state=(filteredF_bW(ii));
    
  if prev_state==1
        prev_m=[1,0];
    else
        prev_m=[0,1];
  end
    
  if outstateF(ii-1)==outstateF(ii-2)
        T=T+1;
    else
        T=1;
        C=C+1;
  end

    markovF=prev_m*(markovF_m^T);
if markovF(1,1)>0.5 
    accA_mk(ii)=markovF(1,1);
    markovF=[1,0];
    p_markov=1;
else
    if markovF(1,2)>0.5
        accF_mk(ii)=markovF(1,2);
        markovF=[0,1];
        p_markov=2;
    else
        p_markov=pred_state;
    end
end

%Decision.

if accF_mk(ii)>=accF_bw(ii)
    outstateF(ii)=p_markov;
else
    outstateF(ii)=pred_state;
end
end

%%
figure(8); 
s(1) = subplot(4, 1, 1);
imagesc(db(A_Spec')); colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Data Fusion "Predicted SW"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*filteredF_sW, 'filled', 'Marker', 'none');
% plot(x, 100.*filtered_bW, 'ob', x, 100.*label_ans, 'sk');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');

s(2) = subplot(4, 1, 2);
imagesc(db(A_Spec')); colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Data Fusion "Predicted BW"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*filteredF_bW, 'filled', 'Marker', 'none');
% plot(x, 100.*filtered_bW, 'ob', x, 100.*label_ans, 'sk');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');


s(3) = subplot(4, 1, 3);
imagesc(db(A_Spec')); colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Data Fusion "Predicted Markov"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*outstateF, 'filled', 'Marker', 'none');
% plot(x, 100.*filtered_bW, 'ob', x, 100.*label_ans, 'sk');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');

s(4) = subplot(4, 1, 4); 
imagesc(db(A_Spec'));% colormap('jet');
hold on
set(gca,'YDir','normal');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 305.5E2 Data Fusion "Correct Label"',...
    'FontSize', 18, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
% stem(x, 100.*label_ans, 'filled', 'Marker', 'none');
% legend({'Predicted Label', 'Correct Label'},...
%     'Location', 'northeast', 'FontSize', 36, 'FontWeight','bold');

%
pF2 = zeros(length(labelsF));

% ===================================== %
% Performance matrix:  p
%
%                   correct label
%                     1        2
%                 -----------------
%              1  - 1/1  -  2/1   -sum
%   predicted     ----------------- 
%              2  - 1/2  -  2/2   -
%                 -----------------
% ===================================== %

for ii = 1:length(label_act)
    if(outstateF(ii) == 1 && label_act(ii) == 1)
        pF2(1,1) = pF2(1,1) +1;
    elseif(outstateF(ii) == 1 && label_act(ii) == 2)
        pF2(1,2) = pF2(1,2) +1;
    elseif(outstateF(ii) == 2 && label_act(ii) == 1)
        pF2(2,1) = pF2(2,1) +1;
    else
        pF2(2,2) = pF2(2,2) +1;
    end
end

pF_2(1,1) = pF2(1,1)/(pF2(1,1)+pF2(2,1));
pF_2(2,1) = pF2(2,1)/(pF2(1,1)+pF2(2,1));
pF_2(1,2) = pF2(1,2)/(pF2(1,2)+pF2(2,2));
pF_2(2,2) = pF2(2,2)/(pF2(1,2)+pF2(2,2));

save('PerformanceMatrixF2.mat', 'pF2', '-mat');


fprintf('After Window Filtering \n');
pF1
fprintf('After Markov Chain \n');
pF2
AccuracyF = (pF2(1,1)+pF2(2,2))/((pF2(1,1)+pF2(1,2)+(pF2(2,1)+pF2(2,2))))