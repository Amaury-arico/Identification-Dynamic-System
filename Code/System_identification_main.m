%% Main - System Identification Model

clear; close all; clc;

%% Load data

fs = 4000; %Hz

% assuming the excited freqeuencies are the same for all realizations and power levels
%excitationpath = fullfile(pwd,'Excitation',filesep);
excitationpath = fullfile(pwd,'Result','2111',strcat('Result_V3_',string(fs),'fs'),filesep);
input_filname = 'Amaury_2111';
load(strcat(excitationpath,input_filname,'_Sel_E0_S0.mat')); % loads variable 'sel'
load(strcat(excitationpath,input_filname,'_Sig_E0_S0.mat')); % loads variable 'td_signal'

[u, y, realizations, power_levels] = acquisition(excitationpath);

N_tot_sample =  size(u,1);                                  % Total nbre sample for 1 realisation
N_sample_per = size(td_signal,1);                           % Nbre sample per periode
Nperiod =  N_tot_sample/N_sample_per;                       %Nbre of periods signal got repeated
freq = (0:N_sample_per-1)*fs/N_sample_per;
excitedbin = sel;

% Differentiate Odd, even and noise for first realisation - NEED FOR ALL
% THE REALISATIONS
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

G = Y./U;

%% Fast Method - 1 Realisation - EXCITED FREQUENCIES ONLY !! - FIND DISTORSIONS LEVELS BY USING U(excited freq at neighbourhood of non excited bins!!)

G_fast = zeros(N_sample_per,1);
G_fast(excitedFreq,1) = Y(excitedFreq,1,1)./U(excitedFreq,1,1);
%Y_odd = Y(oddFreq,1,1);
%Y_even = Y(evenFreq,1,1);
Y_noise = Y(noiseFreq,1,1);

bandinterest = oddFreq<(fs/2);
oddFreq_new = oddFreq(bandinterest==true);
for i=1:length(oddFreq_new)
    if abs(U(oddFreq(i)+6,1,1)) > 30
        new_odd(i) = oddFreq(i)+6;
    elseif abs(U(oddFreq(i)-6,1,1)) > 30
        new_odd(i) = oddFreq(i)-6;
    else
        new_odd(i) = oddFreq(i)-12;
    end
end
Y_odd = Y(oddFreq_new,1,1);
G_fast(oddFreq_new,1) = Y_odd./U(new_odd,1,1);

bandinterest = evenFreq<(fs/2);
evenFreq_new = evenFreq(bandinterest==true);
for i=1:length(evenFreq_new)
    if abs(U(evenFreq(i)+3,1,1)) > 30
        new_even(i) = evenFreq(i)+3;
    elseif abs(U(evenFreq(i)-3,1,1)) > 30
        new_even(i) = evenFreq(i)-3;
    else
        new_even(i) = evenFreq(i)+9;
    end
end
Y_even = Y(evenFreq_new,1,1);
G_fast(evenFreq_new,1) = Y_even./U(new_even,1,1);

bandinterest = noiseFreq<(fs/2);
noiseFreq_new = noiseFreq(bandinterest==true);
for i=1:length(noiseFreq_new)
    if abs(U(noiseFreq(i)+1,1,1)) > 30
        new_noise(i) = noiseFreq(i)+1;
    elseif i > 1 && abs(U(noiseFreq(i)-1,1,1)) > 30
        new_noise(i) = noiseFreq(i)-1;
    elseif abs(U(noiseFreq(i)+2,1,1)) > 30
        new_noise(i) = noiseFreq(i)+2;
    elseif i > 2 && abs(U(noiseFreq(i)-2,1,1)) > 30
        new_noise(i) = noiseFreq(i)-2;
    elseif abs(U(noiseFreq(i)+3,1,1)) > 30
        new_noise(i) = noiseFreq(i)+3;
    elseif abs(U(noiseFreq(i)+4,1,1)) > 30
        new_noise(i) = noiseFreq(i)+4;
    elseif abs(U(noiseFreq(i)+5,1,1)) > 30
        new_noise(i) = noiseFreq(i)+5;
    elseif abs(U(noiseFreq(i)+6,1,1)) > 30
        new_noise(i) = noiseFreq(i)+6;
    elseif abs(U(noiseFreq(i)+7,1,1)) > 30
        new_noise(i) = noiseFreq(i)+7;
    elseif abs(U(noiseFreq(i)+8,1,1)) > 30
        new_noise(i) = noiseFreq(i)+8;
    elseif abs(U(noiseFreq(i)+10,1,1)) > 30
        new_noise(i) = noiseFreq(i)+10;
    elseif abs(U(noiseFreq(i)+11,1,1)) > 30
        new_noise(i) = noiseFreq(i)+11;
    elseif abs(U(noiseFreq(i)+13,1,1)) > 30
        new_noise(i) = noiseFreq(i)+13;
    else
        new_noise(i) = noiseFreq(i)+14;
    end
end
Y_noise = Y(noiseFreq_new,1,1);
G_fast(noiseFreq_new,1) = Y_noise./U(new_noise,1,1);

%% Robust Method -  EXCITED FREQUENCIES ONLY!!! - CAN USE ALL THE FREQUENCIES FOR THE BIN

Nper_effective = Nperiod-1;
G_rea = mean(G((4:6:excitedFreq(end)),:,:),2);
G_rea = squeeze(G_rea);
G_ML = mean(G_rea,2);
for j = 1 : realizations
    add = 0;
    for i = 1 : Nper_effective
        sigma = (G((4:6:excitedFreq(end)),i,j)-G_rea(:,j)).^2 ;
        add = add + sigma;
    end
    sigma_rea(:,j) = add*(1/(Nper_effective*(Nper_effective-1)));
end

mean_variance = (1/realizations^2)*sum(sigma_rea,2);

add = 0;
for j = 1 : realizations
    sigma = (G_rea(:,j)-G_ML).^2 ;
    add = add + sigma;
end
    sigma_ML = add*(1/(realizations*(realizations-1)));




%% Plot FAST Method - 1 realisation

figure;
for k = 1 : realizations
    if mod(realizations,3)==0
    subplot(3,realizations/3,k)
    else
    sublplot(2,realizations/2,k)
    end
    plot(freq(excitedFreq),db(Y(excitedFreq,1,k)),'o', 'LineWidth', 1);
    hold on;
    plot(freq(oddFreq),db(Y(oddFreq,1,k)),'square', 'LineWidth', 1);
    plot(freq(evenFreq),db(Y(evenFreq,1,k)),'+', 'LineWidth', 1);
    plot(freq(noiseFreq),db(Y(noiseFreq,1,k)),'x', 'LineWidth', 1);
    xlim([0 fs/N_sample_per*max(excitedFreq)]);
    hold off;
    xlabel('Frequency [Hz]');
    ylabel('Decibel [dB]');
    legend('Linear','Odd','Even','Noise');
    title(['Single realisation - Realisation n° ',k]);
end

figure;
plot(freq(excitedFreq),db(Y(excitedFreq,1,1)),'o', 'LineWidth', 1);
hold on;
plot(freq(oddFreq),db(Y(oddFreq,1,1)),'square', 'LineWidth', 1);
plot(freq(evenFreq),db(Y(evenFreq,1,1)),'+', 'LineWidth', 1);
plot(freq(noiseFreq),db(Y(noiseFreq,1,1)),'x', 'LineWidth', 1);
xlim([0 fs/N_sample_per*max(excitedFreq)]);
hold off;
xlabel('Frequency [Hz]');
ylabel('Decibel [dB]');
legend('Linear','Odd','Even','Noise');
title('Output Signal - Realisation n° 1');

figure;
subplot(1,2,1)
plot(freq(excitedFreq),db(Y(excitedFreq,1,1)),'o', 'LineWidth', 1);
xlabel('Frequency [Hz]');
ylabel('Decibel [dB]');
title('Single period - Output');
subplot(1,2,2)
plot(freq(excitedFreq),db(U(excitedFreq,1,1)),'o', 'LineWidth', 1);
xlabel('Frequency [Hz]');
ylabel('Decibel [dB]');
title('Single period - Input'); % Reduction of magnitude of the excited signal due to the zero order holder.

figure;

plot(freq(excitedFreq),db(G_fast(excitedFreq,1)),'o', 'LineWidth', 1);
hold on;
plot(freq(oddFreq),db(G_fast(oddFreq,1)),'square', 'LineWidth', 1);
plot(freq(evenFreq),db(G_fast(evenFreq,1)),'+', 'LineWidth', 1);
plot(freq(noiseFreq),db(G_fast(noiseFreq,1)),'x', 'LineWidth', 1);
xlim([0 fs/N_sample_per*max(excitedFreq)]);
title('Fast Method - FRF'); % Reduction of magnitude of the excited signal due to the zero order holder.
xlabel('Frequency [Hz]');
ylabel('Decibel [dB]');
legend('Linear','Odd','Even','Noise');


%% Plot Robust Method

figure;
plot(freq((4:6:excitedFreq(end))),db(G_ML(:,1)),'o', 'LineWidth', 1);
hold on
%plot(freq(oddFreq),db(G_ML(oddFreq,1)),'square', 'LineWidth', 1);
%plot(freq(evenFreq),db(G_ML(evenFreq,1)),'+', 'LineWidth', 1);
%plot(freq(noiseFreq),db(G_ML(noiseFreq,1)),'x', 'LineWidth', 1);
xlim([0 fs/N_sample_per*max(excitedFreq)]);
title('Robust Method')






