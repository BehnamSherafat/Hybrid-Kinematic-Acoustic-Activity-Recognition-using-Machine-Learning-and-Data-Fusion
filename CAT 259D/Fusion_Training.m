
% -----------------------------------------------------------------------
% Thesis
%
% Contact: Behnam Sherafat @ < behnam.sherafat@utah.edu>
% -----------------------------------------------------------------------
%% Initializing
tic;
clear;
time2 = tic;
close all;
clc;

%% Defining STFT Parameters
overlap = 2;
window_size = 2*overlap;
FFTSIZE = 48;
PCA_f = 0.9;


cross_fold = 10;
% delay = 1.7244; % how much audio is prior to vibration
delay = 2.2244;
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

% [A_ori, fsA] = audioread('ZOOM0020.wav');

%% Audio Denoising: method = 'mcra2'; 

% A_ori(A_ori==0)=10^-4;
% A_ori = A_ori(:, 1);
% audiowrite('ZOOM0020_ori.wav', A_ori, fsA);

% Audio Denoising: method = 'mcra2'; 
% specsub_ns('ZOOM0033_ori.wav', 'mcra2', 'ZOOM0033_den.wav')

% Read Modified Audio
[S_Audio, fsA] = audioread('ZOOM0020_den.wav');

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

% Downsample Auido
% down_sample_f = 441;
% S_Audio1 = decimate(S_Audio,441);
% 
% S_Ax1 = S_Ax(1:end-(size(S_Ax,1)-size(S_Audio1,1)));
% S_Ay1 = S_Ay(1:end-(size(S_Ay,1)-size(S_Audio1,1)));
% S_Az1 = S_Az(1:end-(size(S_Az,1)-size(S_Audio1,1)));
% 
% S_Vx1 = S_Vx(1:end-(size(S_Vx,1)-size(S_Audio1,1)));
% S_Vy1 = S_Vy(1:end-(size(S_Vy,1)-size(S_Audio1,1)));
% S_Vz1 = S_Vz(1:end-(size(S_Vz,1)-size(S_Audio1,1)));
% 
% fs1 = fsA/down_sample_f;
% Upsample Vibration
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

%% Plot Figures
% figure(1);
% subplot(2,2,1)
% plot(S_Audio)
% xt1 = get(gca, 'XTick');                                 % 'XTick' Values
% set(gca,'YDir','normal');
% set(gca, 'XTick', xt1, 'XTickLabel', round(xt1/fs2))
% title('CAT 259D Original Audio Signal','FontSize',20,'FontWeight','bold');
% y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
% xlabel('Time (s)','FontSize',18,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
% set(get(gca,'ylabel'),'rotation',0)
% axis tight
% 
% subplot(2,2,3)
% plot(S_Audio1)
% xt2 = get(gca, 'XTick');                                 % 'XTick' Values
% set(gca,'YDir','normal');
% set(gca, 'XTick', xt2, 'XTickLabel', round(xt2/fs1))
% title('CAT 259D Down-sampled Audio Signal','FontSize',20,'FontWeight','bold');
% y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
% xlabel('Time (s)','FontSize',18,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
% set(get(gca,'ylabel'),'rotation',0)
% axis tight
% 
% subplot(2,2,2)
% plot(S_Ax)
% xt3 = get(gca, 'XTick');                                 % 'XTick' Values
% set(gca,'YDir','normal');
% set(gca, 'XTick', xt3, 'XTickLabel', round(xt3/fs1))
% title('CAT 259D Original Acceleration (x) Signal','FontSize',20,'FontWeight','bold');
% y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
% xlabel('Time (s)','FontSize',18,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
% set(get(gca,'ylabel'),'rotation',0)
% axis tight
% 
% subplot(2,2,4)
% plot(S_Ax2)
% xt4 = get(gca, 'XTick');                                 % 'XTick' Values
% set(gca,'YDir','normal');
% set(gca, 'XTick', xt4, 'XTickLabel', round(xt4/fs2))
% title('CAT 259D Up-sampled Acceleration (x) Signal','FontSize',20,'FontWeight','bold');
% y = ylabel('Amplitude','FontSize',18,'FontWeight','bold');
% xlabel('Time (s)','FontSize',18,'FontWeight','bold');
% set(y, 'Units', 'Normalized', 'Position', [-0.15, 0.5, 0]);
% set(get(gca,'ylabel'),'rotation',0)
% axis tight

%% Normalization of Signals

S_Audio = bsxfun(@rdivide,bsxfun(@minus,S_Audio,mean(S_Audio)),std(S_Audio));

S_Ax = bsxfun(@rdivide,bsxfun(@minus,S_Ax,mean(S_Ax)),std(S_Ax));

S_Ay = bsxfun(@rdivide,bsxfun(@minus,S_Ay,mean(S_Ay)),std(S_Ay));

S_Az = bsxfun(@rdivide,bsxfun(@minus,S_Az,mean(S_Az)),std(S_Az));

S_Vx = bsxfun(@rdivide,bsxfun(@minus,S_Vx,mean(S_Vx)),std(S_Vx));

S_Vy = bsxfun(@rdivide,bsxfun(@minus,S_Vy,mean(S_Vy)),std(S_Vy));

S_Vz = bsxfun(@rdivide,bsxfun(@minus,S_Vz,mean(S_Vz)),std(S_Vz));

time3 = toc(time2);

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

short_stop = 1;
arm_raising = 2;
arm_lowering = 2;
shovel_lowering = 2;
shovel_raising = 2;
moving_forward = 3;
moving_backward = 3;
turning_right = 4;
turning_left = 4;

% Adding Labels to Training Data Set

% Add Actual Labels

label_act = zeros(size(Data, 1), 1);
lps = length(label_act)/T; % length per second 

s1 = ceil((0)*lps+1);
f1 = ceil((13.4)*lps);
s2 = ceil((13.4)*lps+1);
f2 = ceil((24.4)*lps);
s3 = ceil((24.4)*lps+1);
f3 = ceil((30.4)*lps);
s4 = ceil((30.4)*lps+1);
f4 = ceil((35.4)*lps);
s5 = ceil((35.4)*lps+1);
f5 = ceil((36.5)*lps);
s6 = ceil((36.5)*lps+1);
f6 = ceil((45.97)*lps);
s7 = ceil((45.97)*lps+1);
f7 = ceil((46.84)*lps);
s8 = ceil((46.84)*lps+1);
f8 = ceil((52.6)*lps);
s9 = ceil((52.6)*lps+1);
f9 = ceil((53.5)*lps);
s10 = ceil((53.5)*lps+1);
f10 = ceil((60.85)*lps);
s11 = ceil((60.85)*lps+1);
f11 = ceil((61.65)*lps);
s12 = ceil((61.65)*lps+1);
f12 = ceil((70.13)*lps);
s13 = ceil((70.13)*lps+1);
f13 = ceil((72.89)*lps);
s14 = ceil((72.89)*lps+1);
f14 = ceil((77.16)*lps);
s15 = ceil((77.16)*lps+1);
f15 = ceil((78.84)*lps);
s16 = ceil((78.84)*lps+1);
f16 = ceil((80.47)*lps);
s17 = ceil((80.47)*lps+1);
f17 = ceil((T)*lps);

label_act(s1:f1) = 2;
label_act(s2:f2) = 2;
label_act(s3:f3) = 2;
label_act(s4:f4) = 2;
label_act(s5:f5) = 2;
label_act(s6:f6) = 3;
label_act(s7:f7) = 1;
label_act(s8:f8) = 3;
label_act(s9:f9) = 1;
label_act(s10:f10) = 4;
label_act(s11:f11) = 1;
label_act(s12:f12) = 4;
label_act(s13:f13) = 2;
label_act(s14:f14) = 4;
label_act(s15:f15) = 2;
label_act(s16:f16) = 2;
label_act(s17:f17) = 1;

% Categorizing Data
data_1 = [Data(s7:f7,:);Data(s9:f9,:);Data(s11:f11,:);Data(s17:f17,:)];
% data_1 = data_1(1:ceil(size(data_1,1)*0.3),:);

data_2 = [Data(s1:f1,:);Data(s2:f2,:);Data(s3:f3,:);Data(s4:f4,:);Data(s5:f5,:);Data(s13:f13,:);Data(s15:f15,:);Data(s16:f16,:)];
% data_2 = data_2(1:ceil(size(data_2,1)*0.3),:);

data_3 = [Data(s6:f6,:);Data(s8:f8,:)];
% data_3 = data_3(1:ceil(size(data_3,1)*0.3),:);

data_4 = [Data(s10:f10,:);Data(s12:f12,:);Data(s14:f14,:)];
% data_4 = data_4(1:ceil(size(data_4,1)*0.3),:);

label_1 = ones(size(data_1,1), 1);
label_2 = 2*ones(size(data_2,1), 1);
label_3 = 3*ones(size(data_3,1), 1);
label_4 = 4*ones(size(data_4,1), 1);

data_set = [data_1;data_2;data_3;data_4];
data_labels = [label_1;label_2;label_3;label_4];

%% 5-fold Cross Validation
% Cross varidation (train: 80%, validate: 10%, test: 10%)
% Shuffle DataSet
for i=1:5
    idx = randperm(size(data_set,1));
    data_set = data_set(idx,:);
    data_labels = data_labels(idx,:);

    cv = cvpartition(size(data_set,1),'HoldOut',0.1);
    idx = cv.test;

    % Test Data
    X_Test  = data_set(idx,:);
    Y_Test  = data_labels(idx,:);

    % Separate labels to train and val data
    X_Train1 = data_set(~idx,:);
    Y_Train1 = data_labels(~idx,:);

    cv = cvpartition(size(X_Train1,1),'HoldOut',0.09);
    idx = cv.test;

    % Validation Data
    X_Val  = X_Train1(idx,:);
    Y_Val  = Y_Train1(idx,:);

    % Train Data
    X_Train = X_Train1(~idx,:);
    Y_Train = Y_Train1(~idx,:);

    model = fitcecoc(X_Train,Y_Train);
    label3 = predict(model,X_Test);

    comparison_fused = table(Y_Test,label3);
    C_fused = confusionmat(Y_Test,label3)
    accuracy_fused = trace(C_fused)/sum(C_fused(:))
    
end

% Cross varidation (train: 90%, validate: 5%, test: 5%)
for i=1:5
    idx = randperm(size(data_set,1));
    data_set = data_set(idx,:);
    data_labels = data_labels(idx,:);

    cv = cvpartition(size(data_set,1),'HoldOut',0.05);
    idx = cv.test;

    % Test Data
    X_Test  = data_set(idx,:);
    Y_Test  = data_labels(idx,:);

    % Separate labels to train and val data
    X_Train1 = data_set(~idx,:);
    Y_Train1 = data_labels(~idx,:);

    cv = cvpartition(size(X_Train1,1),'HoldOut',0.045);
    idx = cv.test;

    % Validation Data
    X_Val  = X_Train1(idx,:);
    Y_Val  = Y_Train1(idx,:);

    % Train Data
    X_Train = X_Train1(~idx,:);
    Y_Train = Y_Train1(~idx,:);

    model = fitcecoc(X_Train,Y_Train);
    label3 = predict(model,X_Test);

    comparison_fused = table(Y_Test,label3);
    C_fused = confusionmat(Y_Test,label3)
    accuracy_fused = trace(C_fused)/sum(C_fused(:))
    
end
%%
% Cross varidation (train: 95%, validate: 2.5%, test: 2.5%)
for i=1:5
    idx = randperm(size(data_set,1));
    data_set = data_set(idx,:);
    data_labels = data_labels(idx,:);

    cv = cvpartition(size(data_set,1),'HoldOut',0.025);
    idx = cv.test;

    % Test Data
    X_Test  = data_set(idx,:);
    Y_Test  = data_labels(idx,:);

    % Separate labels to train and val data
    X_Train1 = data_set(~idx,:);
    Y_Train1 = data_labels(~idx,:);

    cv = cvpartition(size(X_Train1,1),'HoldOut',0.02375);
    idx = cv.test;

    % Validation Data
    X_Val  = X_Train1(idx,:);
    Y_Val  = Y_Train1(idx,:);

    % Train Data
    X_Train = X_Train1(~idx,:);
    Y_Train = Y_Train1(~idx,:);

    model = fitcecoc(X_Train,Y_Train);
    label3 = predict(model,X_Test);

    comparison_fused = table(Y_Test,label3);
    C_fused = confusionmat(Y_Test,label3)
    accuracy_fused = trace(C_fused)/sum(C_fused(:))
    
end

% training_data_1 = Activity_1(1:ceil(size(Activity_1,1)*train_test_ratio),:);
% training_data_2 = Activity_2(1:ceil(size(Activity_2,1)*train_test_ratio),:);
% training_data_3 = Activity_3(1:ceil(size(Activity_3,1)*train_test_ratio),:);
% training_data_4 = Activity_4(1:ceil(size(Activity_4,1)*train_test_ratio),:);
% 
% training_data = [training_data_1;training_data_2;training_data_3;training_data_4];
% 
% testing_data_1 = Activity_1(ceil(size(Activity_1,1)*train_test_ratio)+1:end,:);
% testing_data_2 = Activity_2(ceil(size(Activity_2,1)*train_test_ratio)+1:end,:);
% testing_data_3 = Activity_3(ceil(size(Activity_3,1)*train_test_ratio)+1:end,:);
% testing_data_4 = Activity_4(ceil(size(Activity_4,1)*train_test_ratio)+1:end,:);
% 
% testing_data = [testing_data_1;testing_data_2;testing_data_3;testing_data_4];
% 
% %% Building Fusion Data Test
% 
% fused_training_data = training_data;
% fused_testing_data = testing_data;
% 
% % Implementing PCA
% [COEFF, ~, ~] = pca(fused_training_data(:, 2:end)); % Obtaining coefficients matrix
% p = round(PCA_f*size(fused_training_data, 2)-1);
% Ap = COEFF(:, 1:p);
% 
% training_data_N = [fused_training_data(:, 1), bsxfun(@minus, fused_training_data(:, 2:end), mean(fused_training_data(:, 2:end)))*Ap];
% testing_data_N = [fused_testing_data(:, 1), bsxfun(@minus, fused_testing_data(:, 2:end), mean(fused_training_data(:, 2:end)))*Ap];
% % Creating cross-folds
% % r = randi(cross_fold,size(training_data_N, 1),1);
% % U = [r, training_data_N];
% % firstColumn = U(:, 1);
% % for i = 1:cross_fold
% %      
% %     A = U(firstColumn == i, :);
% %     B = U(firstColumn ~= i, :);
% %     test = A(:, 2:end);
% %     train = B(:, 2:end);
%     
% % SVM Multiclass
% 
% Model = fitcecoc(training_data_N(:,2:end),training_data_N(:,1));
% time9 = tic;
% label3 = predict(Model,testing_data_N(:,2:end));
% time10 = toc(time9)
% comparison_fused = table(testing_data_N(:,1),label3);
% C_fused = confusionmat(testing_data_N(:,1),label3)
% accuracy_fused = trace(C_fused)/sum(C_fused(:))

%% Plotting Confusion Matrices
% 
% figure(8);
% subplot(1,3,1)
% plotconfusion(testing_data_N_1(:,1),label1)
% title('Vibration Confusion Matrix',...
%     'FontSize', 16, 'FontWeight','bold');
% subplot(1,3,2)
% plotconfusion(testing_data_N_2(:,1),label2)
% title('Audio Confusion Matrix',...
%     'FontSize', 16, 'FontWeight','bold');
% subplot(1,3,3)
% plotconfusion(testing_data_N(:,1),label3)
% title('Fused Data Confusion Matrix',...
%     'FontSize', 16, 'FontWeight','bold');
% % toc