function full_sti_signal = generateFullSTISignal(duration, varargin)
% GENERATEFULLSTISIGNAL generates the Full STI test signal for the
% measurement of Speech Transmission Index (STI) according
% to IEC 60268-16 (ed.3) standard using the Full STI method.
%
%   FULL_STI_SIGNAL = GENERATEFULLSTISIGNAL(DURATION) generates a Full STI
%   test signal with each modulation segment lasting DURATION seconds. The
%   total duration includes DURATION seconds for each segment plus silence
%   gaps between segments. Default sampling frequency is 96000 Hz and default
%   silence duration is 0 s.
%
%   FULL_STI_SIGNAL = GENERATEFULLSTISIGNAL(DURATION, FS) generates the Full
%   STI test signal using sampling frequency FS (Hz). The minimum FS is 22050.
%
%   FULL_STI_SIGNAL = GENERATEFULLSTISIGNAL(DURATION, FS, 'silenceDuration',
%   SILENCEDURATION) generates the Full STI test signal with silence gaps 
%   of SILENCEDURATION seconds between each segment (default: 0 s).
%
%   FULL_STI_SIGNAL = GENERATEFULLSTISIGNAL(DURATION, FS, 'doPlot', DOPLOT)
%   when DOPLOT = 1 (default) displays the signal waveform via plotSignalWaveform,
%   and when DOPLOT = 0 suppresses plotting.
%
% References:
%   EN IEC 60268-16:2020 Sound system equipment — Part 16: Objective rating of speech
%   intelligibility by speech transmission index.
%   Záviška, P., Rajmic, P., Schimmel, J. (2024). Matlab Implementation 
%   of STIPA (Speech Transmission Index for Public Address Systems).
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025

% Check number of input arguments
narginchk(1, 8);

    %% Parse input arguments
    [fs, silenceDuration, doPlot] = parseInputs(duration, varargin{:});
    
    %% Generate pink noise
    N = duration * fs;
    pinkNoise = pinknoise(N);
    
    %% Define filtering and modulation settings
    filterOrder = 20;
    octaveBands = [125, 250, 500, 1000, 2000, 4000, 8000];                % Octave band center frequencies
    fm = [0.63, 0.8, 1, 1.25, 1.6, 2, 2.5, 3.15, 4, 5, 6.3, 8, 10, 12.5]; % Modulation frequencies (Hz)
    m = 1;               % Full modulation depth
    t = (0:N-1).' / fs;  % Time vector
    
    % Speech spectrum levels (dB) per band, converted to linear scale
    levels = [-2.5, 0.5, 0, -6, -12, -18, -24];
    G = 10 .^ (levels / 20);
    
    %% Generate filtered pink noise for each octave band
    filteredPinkNoise = generateFilteredPinkNoise(pinkNoise, fs, octaveBands, filterOrder);
    
    %% Generate modulated signals and assemble the final signal
    [full_sti_signal, activeRMS] = generateModulatedSignals(filteredPinkNoise, fm, m, t, N, fs, octaveBands, G, silenceDuration);
   
    % Normalization for the Full STI signal
    targetRMS     = 0.07;                    % Target RMS level for the active signal
    scalingFactor = targetRMS / activeRMS;   % Calculate the scaling factor

    % Apply scaling to the entire vector — zeros remain zeros
    full_sti_signal = full_sti_signal * scalingFactor;

    %% Plot signal waveform if requested
    if doPlot
        titleStr = sprintf('Full STI: %ds, gap %ds at %d Hz', duration, silenceDuration, fs);
        plotSignalWaveform(full_sti_signal, fs, titleStr);
    end
end

%% Helper Functions

function [fs, silenceDuration, doPlot] = parseInputs(duration, varargin)
% PARSEINPUTS Parses input arguments and returns sampling frequency,
% silence duration, and plot option.
    p = inputParser;
    validScalarPosNum    = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    validSamplingFreq    = @(x) isnumeric(x) && isscalar(x) && (x >= 22050);
    validSilenceDuration = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    validPlotOption      = @(x) isnumeric(x) && isscalar(x) && (x == 0 || x == 1);
    
    addRequired(p, 'duration', validScalarPosNum);
    
    defaultFs = 96000;
    addOptional(p, 'fs', defaultFs, validSamplingFreq);
    
    defaultSilenceDuration = 0;
    addOptional(p, 'silenceDuration', defaultSilenceDuration, validSilenceDuration);
    
    defaultDoPlot = 1;
    addOptional(p, 'doPlot', defaultDoPlot, validPlotOption);
    
    parse(p, duration, varargin{:});
    fs     = p.Results.fs;
    doPlot = p.Results.doPlot;
    silenceDuration = p.Results.silenceDuration;
end

function filteredPinkNoise = generateFilteredPinkNoise(pinkNoise, fs, octaveBands, filterOrder)
%GENERATEFILTEREDPINKNOISE Filters the input pink noise for each octave band.
    N = length(pinkNoise);
    numBands = length(octaveBands);
    filteredPinkNoise = zeros(N, numBands);
    for idx = 1:numBands
        octFilt = octaveFilter(octaveBands(idx), '1/2 octave', ...
            'SampleRate', fs, 'FilterOrder', filterOrder);
        filteredPinkNoise(:, idx) = octFilt(pinkNoise);
    end
end

function [signal, activeRMS] = generateModulatedSignals(filteredPinkNoise, fm, m, t, N, fs, octaveBands, G, silenceDuration)
%GENERATEMODULATEDSIGNALS Generates modulated signals for each modulation
%frequency and octave band, and assembles them into the final signal.
    numBands = length(octaveBands);
    numModFreq = length(fm);
    totalSignals = numBands * numModFreq;
    
    silence_gap = round(silenceDuration * fs);  % Number of samples for silence gap
    totalNumSamples = totalSignals * N + (totalSignals - 1) * silence_gap;
    signal = zeros(totalNumSamples, 1);
    
    currentIndex = 1;
    % signalCount = 1;
    sumSq = 0;
    for fmIdx = 1:numModFreq
        % Compute the modulation signal once for the current modulation frequency
        modVec = sqrt(0.5 * (1 + m * cos(2 * pi * fm(fmIdx) * t)));
        for bandIdx = 1:numBands            
            % Apply modulation to the pre-filtered pink noise for the current band
            modulatedSignal = filteredPinkNoise(:, bandIdx) .* modVec;
            modulatedSignal = modulatedSignal * G(bandIdx);  % Apply band level
            sumSq = sumSq + sum(modulatedSignal.^2);
            
            % Insert the modulated signal into the final signal with silence gap
            signal(currentIndex:currentIndex+N-1) = modulatedSignal;
            currentIndex = currentIndex + N + silence_gap;
            % signalCount = signalCount + 1;
        end
    end
    activeRMS = sqrt(sumSq / (N * totalSignals));
end