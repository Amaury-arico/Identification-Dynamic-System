%% Result Export

function [u, y, realizations, power_levels] = acquisition(path)
    %filename = strcat('');
    [realizations, power_levels] = detect_counts(path);
    
    % load a single file to get N
    file = strcat(path, "ACQ_R0_P0_E0_M0_F0.mat");
    disp(file);
    load(file);
    N = length(YR0);

    for r = 0:realizations-1
        for p = 0:power_levels-1
            file = strcat(path, "ACQ_R", num2str(r), "_P", num2str(p), "_E0_M0_F0.mat");
            load(file);
            u(:, r+1, p+1) = YR0(:);
            y(:, r+1, p+1) = YR1(:);
        end
    end
end

function [realizations, power_levels] = detect_counts(resultsPath)
% Detect number of realizations and power levels by probing for existing files.
% Returns counts (non-negative integers). Filenames checked:
% [resultsPath]/[filename]ACQ_R{r}_P{p}_E0_M0_F0.mat
%
% Example: [realizations, power_levels] = detect_counts('../results/', 'out2k');
% detect realizations by probing P0 files starting at R0
    r = 0;
    while true
        fname = fullfile(resultsPath, strcat("ACQ_R", num2str(r), "_P0_E0_M0_F0.mat"));
        if isfile(fname)
            r = r + 1;
        else
            break;
        end
    end
    realizations = r; % may be 0 if no files found

    % detect power levels by probing R0 files starting at P0 (only if at least one realization)
    p = 0;
    if realizations > 0
        while true
            fname = fullfile(resultsPath, strcat("ACQ_R0_P", num2str(p), "_E0_M0_F0.mat"));
            if isfile(fname)
                p = p + 1;
            else
                break;
            end
        end
    end
    power_levels = p; % may be 0 if no files found
end