%% Main - System Identification Model

clear; close all; clc;

%% Load data

fs = 1000; %Hz

% assuming the excited freqeuencies are the same for all realizations and power levels
%excitationpath = fullfile(pwd,'Excitation',filesep);
excitationpath = fullfile(pwd,'Result','2111','Result_V2_1000fs',filesep);
input_filname = 'Amaury_2111';
load(strcat(excitationpath,input_filname,'_Sel_E0_S0.mat')); % loads variable 'sel'
load(strcat(excitationpath,input_filname,'_Sig_E0_S0.mat')); % loads variable 'td_signal'

[u, y, realizations, power_levels] = acquisition(excitationpath);

N_tot_sample =  size(u,1); % Total nbre sample for 1 realisation
N_sample_per = size(td_signal,1); % Nbre sample per periode
Nperiod =  N_tot_sample/N_sample_per; %Nbre of periods signal got repeated
freq = (0:N_sample_per-1)*fs/N_sample_per;
excitedbin = sel;

% Differentiate Odd, even and noise for first realisation
excitedFreq = (excitedbin(:, 1) + 1)';
oddFreq = setdiff((4:6:N_sample_per), excitedFreq);
evenFreq = setdiff((7:6:(N_sample_per)), excitedFreq);
noiseFreq = setdiff(1:N_sample_per, [excitedFreq, oddFreq, evenFreq]);

% Remove 1 period of measurement - Transient
Transient_period = 1;
u_period = reshape(u, N_sample_per,Nperiod,realizations);
u_period_no_transient = u_period(:,1+Transient_period:end,:);
y_period = reshape(y, N_sample_per,Nperiod,realizations);
y_period_no_transient = y_period(:,1+Transient_period:end,:);

for r = 1 : realizations
    for i = 1 : Nperiod-Transient_period
        U(:,i,r) = fft(u_period_no_transient(:,i,r));
        Y(:,i,r) = fft(y_period_no_transient(:,i,r));
    end
end

%G = mean(Y(:,:,1),2)./mean(U(:,:,1),2);
%G_per = G;
G = Y./U;
G_per = mean(G(:,:,1),2);

%% Plot 1 period - realisation

figure;
plot(freq(excitedFreq),db(G_per(excitedFreq,1)),'o', 'LineWidth', 1);
hold on
plot(freq(oddFreq),db(G_per(oddFreq,1)),'square', 'LineWidth', 1);
%plot(freq(evenFreq),db(G_per(evenFreq,1)),'+', 'LineWidth', 1);
%plot(freq(noiseFreq),db(G_per(noiseFreq,1)),'x', 'LineWidth', 1);
xlim([0 fs/N_sample_per*max(excitedFreq)]);
title('Single realisation')







