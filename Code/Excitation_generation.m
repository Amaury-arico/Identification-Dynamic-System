%% Identification of dynamical system - Excitation generation

% Name : Amaury Arico
% University - Department : Bruface ULB-VUB
% Option : MA2-IRELE - Telecom/control
% Date : 20251114

% System : Silver box

%Login for the PXI computers:
%User = .\student
%Password = Student

clear; close all;
%% Excitation

fs = 3000;                       % Hz
N = 3000;                        % Number of samples (nbre of bins)
freq = (0:N-1)*fs/N;             % Frequencies axis
t = (0:N-1)/fs;                  % Time axis

M = 12;                          % Realisation
Q = 6;                          % Signal repeat
A = 1;                           % Magnitude (V)

rms_desired = 0.1;

fft_signal = zeros(N,M);             % Signal matrix
td_signal = zeros(N,M);

%% Excitation of N/2 bins (0 to 750 Hz)

Odd_start_bin = 4;  % Bin n° 3 - let 2 bin for noise captures
it = (Odd_start_bin-1)*2; % Spearation between 2 odd bins
random_hole = 0;
r=0;
residmaxbin = mod((N/4)+1,(3*it)+4);
maxbin = (N/4)+1+residmaxbin;
nbrebinmax = round(((N/4)-3)/it);
sel = zeros(nbrebinmax-round(maxbin/(3*it))-1,M);
j = 1;
hole = zeros(M,1);

for m = 1:M % 12 realisation with random phases - NL behave like a noisy disturbances
    j = 1;
    r=0;
    for k = Odd_start_bin:it:maxbin

        % Select bin for holes
        if mod(k-4,(3*it)) == 0 
            random_hole = round(2*rand);
            r=r+1;
        end

        % Add holes
        if mod(k-4,(3*it)) ~= random_hole*it
            phi = unifrnd(-pi,pi);             % Random phase
            fft_signal(k,m) = A*exp(1j*phi);   % TD = 1/N sum(FTT(k)*exp(-j*2pi*k*n/N)
            sel(j,m) = k-1;
            j=j+1;
        else
            hole(m,1) = hole(m,1) +1;
        end
    end

    td_signal(:,m) = real(ifft(fft_signal(:,m)));
    td_signal(:,m) = td_signal(:,m)*rms_desired/rms(td_signal(:,m));
    %full_td_signal(:,(1+(m-1)*Q):m*Q) = repmat(td_signal(:,m),1,Q); 
end

full_td_signal(:, :) = repmat(td_signal(:,:), Q, 1);

fd_signal = fft(td_signal);
full_fd_signal = fft(full_td_signal);

resultsPath = fullfile(pwd,'Excitation',filesep); 
filename = 'Amaury_2111';

save([resultsPath,filename '_Sig_E0_S0.mat'], 'full_td_signal');
save([resultsPath,filename '_Sel_E0_S0.mat'], 'sel');

figure;
for k = 1:6
    subplot(2,3,k);
    plot(freq(1:42), db(fd_signal(1:42,k)));
    hold on;
    title('Excitation - 6 realisations');
end
hold off;









