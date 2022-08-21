% -----------------------------------------------------------------------
% Thesis
%
% Contact: Behnam Sherafat @ < behnam.sherafat@utah.edu>
% -----------------------------------------------------------------------
%% Initializing
tic
clear;
close all;
clc;

%% Defining STFT Parameters
overlap = 256;
window_size = 2*overlap;
FFTSIZE = 48;
PCA_f = 0.9;
cross_fold = 10;
delay = 0.07; % how much audio is prior to vibration
train_test_ratio = 0.9;

%% Reading Sensor Data

load test2.mat

fsAccel = length(timestampAccel)/timestampAccel(length(timestampAccel));
fsAng = length(timestampAng)/timestampAng(length(timestampAng));

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

S_Ax = Refilled_Ax;
S_Ay = Refilled_Ay;
S_Az = Refilled_Az;

% Read other Sensor Data
S_Vx = Vx;
S_Vy = Vy;
S_Vz = Vz;

%% Reading Audio Data

[A_ori, fsA] = audioread('ZOOM0054.wav');

%% Audio Denoising: method = 'mcra2'; 

A_ori(A_ori==0)=10^-4;
A_ori = A_ori(:, 1);
audiowrite('ZOOM0054_ori.wav', A_ori, fsA);

% Audio Denoising: method = 'mcra2'; 
% specsub_ns('ZOOM0054_ori.wav', 'mcra2', 'ZOOM0054_den.wav')

% Read Modified Audio
[S_Audio, fsA] = audioread('ZOOM0054_den.wav');

%% Mutual Time Stamp and Sycnhronization
start_Audio = round(delay*fsA);
S_Audio = S_Audio(start_Audio:end);
start_Audio = 0;
T_A = size(S_Audio, 1)/fsA;
% 
T_Accel = size(S_Ax, 1)/fsAccel;
T_Ang = size(S_Vx, 1)/fsAng;

%% Resampling Vibration Signals

T1 = [T_Accel T_Ang T_A, 120];
T = floor(min(T1)); % in seconds

S_Ax = S_Ax(1:floor(T*fsAccel)+1);
S_Ay = S_Ay(1:floor(T*fsAccel)+1);
S_Az = S_Az(1:floor(T*fsAccel)+1);

S_Vx = S_Vx(1:floor(T*fsAng)+1);
S_Vy = S_Vy(1:floor(T*fsAng)+1);
S_Vz = S_Vz(1:floor(T*fsAng)+1);

S_Audio = S_Audio(1:floor(T*fsA)+1);

xA = linspace(0, T, size(S_Ax,1));
xV = linspace(0, T, size(S_Vx,1));
xq = 0:1/fsA:T;

S_Ax = (interp1(xA, S_Ax, xq, 'linear'))';
S_Ay = (interp1(xA, S_Ay, xq, 'linear'))';
S_Az = (interp1(xA, S_Az, xq, 'linear'))';

S_Vx = (interp1(xV, S_Vx, xq, 'linear'))';
S_Vy = (interp1(xV, S_Vy, xq, 'linear'))';
S_Vz = (interp1(xV, S_Vz, xq, 'linear'))';

fs = fsA;

%% Normalization of Signals

S_Audio = bsxfun(@rdivide,bsxfun(@minus,S_Audio,mean(S_Audio)),std(S_Audio));

S_Ax = bsxfun(@rdivide,bsxfun(@minus,S_Ax,mean(S_Ax)),std(S_Ax));

S_Ay = bsxfun(@rdivide,bsxfun(@minus,S_Ay,mean(S_Ay)),std(S_Ay));

S_Az = bsxfun(@rdivide,bsxfun(@minus,S_Az,mean(S_Az)),std(S_Az));

S_Vx = bsxfun(@rdivide,bsxfun(@minus,S_Vx,mean(S_Vx)),std(S_Vx));

S_Vy = bsxfun(@rdivide,bsxfun(@minus,S_Vy,mean(S_Vy)),std(S_Vy));

S_Vz = bsxfun(@rdivide,bsxfun(@minus,S_Vz,mean(S_Vz)),std(S_Vz));

%% Building Vibration and Audio Data Set

% ============================== Training Ax ============================== %
[V_Spec1, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Ax, 0, T, fs, FFTSIZE, overlap, window_size);
F_Ax = [V_Spec1, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

% ============================== Training Ay ============================== %
[V_Spec2, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Ay, 0, T, fs, FFTSIZE, overlap, window_size);
F_Ay = [V_Spec2, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

% ============================== Training Az ============================== %
[V_Spec3, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Az, 0, T, fs, FFTSIZE, overlap, window_size);
F_Az = [V_Spec3, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

% ============================== Training Vx ============================== %
[V_Spec4, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vx, 0, T, fs, FFTSIZE, overlap, window_size);
F_Vx = [V_Spec4, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

% ============================== Training Vy ============================== %
[V_Spec5, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vy, 0, T, fs, FFTSIZE, overlap, window_size);
F_Vy = [V_Spec5, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

% ============================== Training Vz ============================== %
[V_Spec6, V_SpecM, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO] = FunctionExtractionV(S_Vz, 0, T, fs, FFTSIZE, overlap, window_size);
F_Vz = [V_Spec6, V_RMS, V_E, V_ZCR, V_SF, V_SpecEntropy, V_SC, V_SRO];

% ============================== Training Audio ============================== %
[A_Spec, A_SpecM, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO] = FunctionExtractionA(S_Audio, 0, T, fs, FFTSIZE, overlap, window_size);
F_A = [A_Spec, A_RMS, A_E, A_SF, A_SpecEntropy, A_SC, A_SRO];

% ==================================================================== %

Data = [abs(F_Ax) abs(F_Ay) abs(F_Az) abs(F_Vx) abs(F_Vy) abs(F_Vz) abs(F_A)]; % concatenate all STFT matrix
Data = bsxfun(@rdivide,bsxfun(@minus,Data,mean(Data)),std(Data));

%% Defining Classes
%dozer staying still with blade lowered 1
%moving forward with blade lowered 2
%signal going off 1
%quick stop 1
%moving forward with blade lowered 2
%quick stop and blade raised 1
%moving backwards 3
%quick stop and blade lowered 1
%moving forward with blade lowered 2

% Adding Labels to Training Data Set

% Add Actual Labels

label_act = zeros(size(Data, 1), 1);
lps = length(label_act)/T; % length per second 

s1 = ceil((0)*lps+1);
f1 = ceil((1)*lps);
s2 = ceil((1)*lps+1);
f2 = ceil((12)*lps);
s3 = ceil((12)*lps+1);
f3 = ceil((14)*lps);
s4 = ceil((14)*lps+1);
f4 = ceil((50)*lps);
s5 = ceil((50)*lps+1);
f5 = ceil((51)*lps);
s6 = ceil((51)*lps+1);
f6 = ceil((89)*lps);
s7 = ceil((89)*lps+1);
f7 = ceil((90)*lps);
s8 = ceil((90)*lps+1);
f8 = ceil((T)*lps);

label_act(s1:f1) = 1;
label_act(s2:f2) = 2;
label_act(s3:f3) = 1;
label_act(s4:f4) = 2;
label_act(s5:f5) = 1;
label_act(s6:f6) = 3;
label_act(s7:f7) = 1;
label_act(s8:f8) = 2;

% Categorizing Data
Activity_1 = [Data(s7:f7,:);Data(s3:f3,:);Data(s5:f5,:);Data(s7:f7,:)];
Activity_2 = [Data(s2:f2,:);Data(s4:f4,:);Data(s8:f8,:)];
Activity_3 = Data(s6:f6,:);

X = [Activity_1; Activity_2; Activity_3];
Y = [ones(size(Activity_1,1), 1); 2*ones(size(Activity_2,1), 1); 3*ones(size(Activity_3,1), 1)];

[idx,scores] = fscchi2(X,Y);
idxInf = find(isinf(scores));
find(isinf(scores))
bar(scores(idx))
xlabel('Predictor rank')
ylabel('Predictor importance score')
hold on
bar(scores(idx(length(idxInf)+1))*ones(length(idxInf),1))
legend('Finite Scores','Inf Scores')
hold off
set(gca,'fontsize', 16) 


Activity_1 = [ones(size(Activity_1,1), 1), Activity_1];
Activity_2 = [2*ones(size(Activity_2,1), 1), Activity_2];
Activity_3 = [3*ones(size(Activity_3,1), 1), Activity_3];

Activity_1 = Activity_1(randperm(size(Activity_1,1)),:);
Activity_2 = Activity_2(randperm(size(Activity_2,1)),:);
Activity_3 = Activity_3(randperm(size(Activity_3,1)),:);

training_data_1 = Activity_1(1:ceil(size(Activity_1,1)*train_test_ratio),:);
training_data_2 = Activity_2(1:ceil(size(Activity_2,1)*train_test_ratio),:);
training_data_3 = Activity_3(1:ceil(size(Activity_3,1)*train_test_ratio),:);

training_data = [training_data_1;training_data_2;training_data_3];

testing_data_1 = Activity_1(ceil(size(Activity_1,1)*train_test_ratio)+1:end,:);
testing_data_2 = Activity_2(ceil(size(Activity_2,1)*train_test_ratio)+1:end,:);
testing_data_3 = Activity_3(ceil(size(Activity_3,1)*train_test_ratio)+1:end,:);

testing_data = [testing_data_1;testing_data_2;testing_data_3];

%% Plotting Data

x = 1:length(label_act);

figure(1);
subplot(2,1,1)
plot(S_Ax)
xt1 = get(gca, 'XTick');                                 % 'XTick' Values
set(gca,'YDir','normal');
set(gca, 'XTick', xt1, 'XTickLabel', round(xt1/fs))
title('CAT 259D Ax Signal','FontSize',20,'FontWeight','bold');
y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight
subplot(2,1,2)
imagesc(db(V_Spec1')); colormap('jet'); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 259D Ax Signal Spectrogram',...
    'FontSize', 20, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
axis tight
hold on
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
axis tight

figure(2);
subplot(2,1,1)
plot(S_Ay)
xt1 = get(gca, 'XTick');                                 % 'XTick' Values
set(gca,'YDir','normal');
set(gca, 'XTick', xt1, 'XTickLabel', round(xt1/fs))
title('CAT 259D Ay Signal','FontSize',20,'FontWeight','bold');
y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight
subplot(2,1,2)
imagesc(db(V_Spec2')); colormap('jet'); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 259D Ay Signal Spectrogram',...
    'FontSize', 20, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
axis tight
hold on
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
axis tight


figure(3);
subplot(2,1,1)
plot(S_Az)
xt1 = get(gca, 'XTick');                                 % 'XTick' Values
set(gca,'YDir','normal');
set(gca, 'XTick', xt1, 'XTickLabel', round(xt1/fs))
title('CAT 259D Az Signal','FontSize',20,'FontWeight','bold');
y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight
subplot(2,1,2)
imagesc(db(V_Spec3')); colormap('jet'); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 259D Az Signal Spectrogram',...
    'FontSize', 20, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
axis tight
hold on
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
axis tight


figure(4);
subplot(2,1,1)
plot(S_Vx)
xt1 = get(gca, 'XTick');                                 % 'XTick' Values
set(gca,'YDir','normal');
set(gca, 'XTick', xt1, 'XTickLabel', round(xt1/fs))
title('CAT 259D Vx Signal','FontSize',20,'FontWeight','bold');
y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight
subplot(2,1,2)
imagesc(db(V_Spec4')); colormap('jet'); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 259D Vx Signal Spectrogram',...
    'FontSize', 20, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
axis tight
hold on
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
axis tight


figure(5);
subplot(2,1,1)
plot(S_Vy)
xt1 = get(gca, 'XTick');                                 % 'XTick' Values
set(gca,'YDir','normal');
set(gca, 'XTick', xt1, 'XTickLabel', round(xt1/fs))
title('CAT 259D Vy Signal','FontSize',20,'FontWeight','bold');
y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight
subplot(2,1,2)
imagesc(db(V_Spec5')); colormap('jet'); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 259D Vy Signal Spectrogram',...
    'FontSize', 20, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
axis tight
hold on
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
axis tight


figure(6);
subplot(2,1,1)
plot(S_Vz)
xt1 = get(gca, 'XTick');                                 % 'XTick' Values
set(gca,'YDir','normal');
set(gca, 'XTick', xt1, 'XTickLabel', round(xt1/fs))
title('CAT 259D Vz Signal','FontSize',20,'FontWeight','bold');
y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight
subplot(2,1,2)
imagesc(db(V_Spec6')); colormap('jet'); set(gca,'YDir','normal'); c = get(gca,'Clim');
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
title('CAT 259D Vz Signal Spectrogram',...
    'FontSize', 20, 'FontWeight','bold');
y=ylabel({'Normalized', 'Frequency'}, 'FontSize', 18, 'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
xlabel('Time (s)', 'FontSize', 10, 'FontWeight','bold');
grid off
axis tight
hold on
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
axis tight


figure(7);
subplot(2,1,1)
plot(S_Audio)
xt2 = get(gca, 'XTick');                                 % 'XTick' Values
set(gca,'YDir','normal');
set(gca, 'XTick', xt2, 'XTickLabel', round(xt2/fs))
title('CAT 259D Audio Signal','FontSize',20,'FontWeight','bold');
y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
xlabel('Time (s)','FontSize',10,'FontWeight','bold');
set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(get(gca,'ylabel'),'rotation',0)
axis tight
subplot(2,1,2)
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
hold on
stem(x, (FFTSIZE/8).*label_act, 'filled', 'Marker', 'none', 'Color', 'black');
axis tight


%% Building Vibration Data Fusion Test

vibration_training_data = training_data(:,1:((((FFTSIZE/2)+1)+7)*6)+1);
vibration_testing_data = testing_data(:,1:((((FFTSIZE/2)+1)+7)*6)+1);

% Implementing PCA
[COEFF, ~, ~] = pca(vibration_training_data(:, 2:end)); % Obtaining coefficients matrix
p = round(PCA_f*size(vibration_training_data, 2)-1);
Ap = COEFF(:, 1:p);

training_data_N_1 = [vibration_training_data(:, 1), bsxfun(@minus, vibration_training_data(:, 2:end), mean(vibration_training_data(:, 2:end)))*Ap];
testing_data_N_1 = [vibration_testing_data(:, 1), bsxfun(@minus, vibration_testing_data(:, 2:end), mean(vibration_training_data(:, 2:end)))*Ap];
% Creating cross-folds
% r = randi(cross_fold,size(training_data_N, 1),1);
% U = [r, training_data_N];
% firstColumn = U(:, 1);
% for i = 1:cross_fold
%      
%     A = U(firstColumn == i, :);
%     B = U(firstColumn ~= i, :);
%     test = A(:, 2:end);
%     train = B(:, 2:end);
    
% SVM Multiclass

Model = fitcecoc(training_data_N_1(:,2:end),training_data_N_1(:,1));
label = predict(Model,testing_data_N_1(:,2:end));
comparison_vibration = table(testing_data_N_1(:,1),label);
C_vibration = confusionmat(testing_data_N_1(:,1),label)
accuracy_vibration = trace(C_vibration)/sum(C_vibration(:))

%% Building Audio Data Test

audio_training_data = [training_data(:,1), training_data(:,((((FFTSIZE/2)+1)+7)*6)+2:end)];
audio_testing_data = [testing_data(:,1), testing_data(:,((((FFTSIZE/2)+1)+7)*6)+2:end)];

% Implementing PCA
[COEFF, ~, ~] = pca(audio_training_data(:, 2:end)); % Obtaining coefficients matrix
p = round(PCA_f*size(audio_training_data, 2)-1);
Ap = COEFF(:, 1:p);

training_data_N_2 = [audio_training_data(:, 1), bsxfun(@minus, audio_training_data(:, 2:end), mean(audio_training_data(:, 2:end)))*Ap];
testing_data_N_2 = [audio_testing_data(:, 1), bsxfun(@minus, audio_testing_data(:, 2:end), mean(audio_training_data(:, 2:end)))*Ap];
% Creating cross-folds
% r = randi(cross_fold,size(training_data_N, 1),1);
% U = [r, training_data_N];
% firstColumn = U(:, 1);
% for i = 1:cross_fold
%      
%     A = U(firstColumn == i, :);
%     B = U(firstColumn ~= i, :);
%     test = A(:, 2:end);
%     train = B(:, 2:end);
    
% SVM Multiclass

Model = fitcecoc(training_data_N_2(:,2:end),training_data_N_2(:,1));
label = predict(Model,testing_data_N_2(:,2:end));
comparison_audio = table(testing_data_N_2(:,1),label);
C_audio = confusionmat(testing_data_N_2(:,1),label)
accuracy_audio = trace(C_audio)/sum(C_audio(:))

%% Building Fusion Data Test

fused_training_data = training_data;
fused_testing_data = testing_data;

% Implementing PCA
[COEFF, ~, ~] = pca(fused_training_data(:, 2:end)); % Obtaining coefficients matrix
p = round(PCA_f*size(fused_training_data, 2)-1);
Ap = COEFF(:, 1:p);

training_data_N = [fused_training_data(:, 1), bsxfun(@minus, fused_training_data(:, 2:end), mean(fused_training_data(:, 2:end)))*Ap];
testing_data_N = [fused_testing_data(:, 1), bsxfun(@minus, fused_testing_data(:, 2:end), mean(fused_training_data(:, 2:end)))*Ap];
% Creating cross-folds
% r = randi(cross_fold,size(training_data_N, 1),1);
% U = [r, training_data_N];
% firstColumn = U(:, 1);
% for i = 1:cross_fold
%      
%     A = U(firstColumn == i, :);
%     B = U(firstColumn ~= i, :);
%     test = A(:, 2:end);
%     train = B(:, 2:end);
    
% SVM Multiclass

Model = fitcecoc(training_data_N(:,2:end),training_data_N(:,1));
label = predict(Model,testing_data_N(:,2:end));
comparison_fused = table(testing_data_N(:,1),label);
C_fused = confusionmat(testing_data_N(:,1),label)
accuracy_fused = trace(C_fused)/sum(C_fused(:))

toc