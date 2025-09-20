 function [STI, mk] = stipa(signal, fs, varargin)
% STIPA calculates the Speech Transmission Index (STI) from a recorded 
% STIPA test signal. The calculation follows the IEC 60268-16 (ed.3) 
% standard using the direct STIPA % (Speech Transmission Index for Public 
% Address Systems) method.
%
%   [STI, MK] = STIPA(SIGNAL, FS) computes the Speech Transmission Index
%   from the test signal SIGNAL and its sampling frequency FS in Hz,
%   returning the Speech Transmission index STI and a 2-by-7 matrix of the
%   respective modulation transfer values MK for each octave band and
%   modulation frequency.
%
%   [STI, MK] = STIPA(SIGNAL, FS, REFERENCE) computes the Speech
%   Transmission Index using the REFERENCE signal instead of the default
%   value 0.55, which results from the STIPA test signal. If the sampling
%   frequency of the REFERENCE is not specified, it is assumed that the
%   sampling frequency is the same as the sampling frequency of the test
%   SIGNAL FS.
%
%   [STI, MK] = STIPA(SIGNAL, FS, REFERENCE, FSREF) computes the Speech
%   Transmission Index using the REFERENCE signal and its sampling
%   frequency FSREF.
%   
%   [STI, MK] = STIPA(SIGNAL, FS, 'Lsk', LSK) computes the Speech
%   Transmission Index with adjustment of the MTF for auditory masking and 
%   threshold effects.
%
%   [STI, MK] = STIPA(SIGNAL, FS, 'Lsk', LSK, 'Lnk', LNK) computes the
%   Speech Transmission Index with adjustments of the MTF for ambient 
%   noise, and auditory masking and threshold effects.
%
%   [STI, MK] = STIPA(SIGNAL, FS, 'SignalStart', SIGNALSTART) manually
%   specifies the start of the test signal in the recording (in seconds).
%   If not provided, the start is automatically detected using an energy
%   threshold.
%
%   [STI, MK] = STIPA(SIGNAL, FS, 'doPlot', DOPLOT) when 'doPlot' is set to 1
%   (default), the function displays graphical output of the results.
%   In this mode, the user is prompted to select an optional calibration
%   reference signal (94 or 114 dB SPL) to compute absolute dB SPL levels.
%   If the prompt is canceled, levels are shown relative to the input units.
%   When doPlot is set to 0, the graphical output is suppressed.
%
%   [STI, MK] = STIPA(SIGNAL, FS, 'doTable', DOTABLE) displays 
%   the MTF matrix in a table when DOTABLE is set to 1 (default). 
%   If DOTABLE is set to 0, the table is not shown.
%
% Reference:
%   EN IEC 60268-16:2020 Sound system equipment — Part 16: Objective rating of speech
%   intelligibility by speech transmission index.
%
% Copyright Pavel Záviška, Brno University of Technology, 2023-2024
% Based on original code by Pavel Záviška; modified version by Šimon Cieslar

    % Check number of input arguments
    narginchk(2, 14);
    
    % Parse input arguments:
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    validColumnVector = @(x) isnumeric(x) && iscolumn(x);
    valid7PosVector   = @(x) isnumeric(x) && length(x) == 7 && all(x > 0);
    validBinaryFlag   = @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1);
    validNonNegativeScalar = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    
    addRequired(p, 'signal', validColumnVector);
    addRequired(p, 'fs', validScalarPosNum);
    addOptional(p, 'reference', NaN, validColumnVector);
    addOptional(p, 'fsRef', fs, validScalarPosNum);
    addParameter(p, 'Lsk', NaN, valid7PosVector);
    addParameter(p, 'Lnk', NaN, valid7PosVector);
    addParameter(p, 'doPlot', 1, validBinaryFlag);
    addParameter(p, 'doTable', 1, validBinaryFlag);
    addParameter(p, 'SignalStart', NaN, validNonNegativeScalar);
    parse(p, signal, fs, varargin{:});

    % Definition of weighting and redundancy factors
    % weighting factors for male speech according to the standard
    % [ page 48 - Annex A ]
    alpha_k = [0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173];
    
    % redundancy factors for male speech according to the standard
    % [ page 48 - Annex A ]
    beta_k = [0.085, 0.078, 0.065, 0.011, 0.047, 0.095];
        
    % Determine the start of the signal (manual vs. automatic)
    if isnan(p.Results.SignalStart)
        % Automatic detection using an energy threshold
        defaultThreshold = 0.01 * max(signal.^2); % threshold can be adjusted as needed
        startIndex = detectSignalStart(signal, fs, defaultThreshold);
        fprintf('Automatic detection: Signal starts approximately at %.2f s.\n', startIndex/fs);
    else
        startIndex = p.Results.SignalStart * fs;
        if startIndex < 1
            startIndex = 1;
        end
        fprintf('Manually provided signal start: %.2f s.\n', p.Results.SignalStart);
    end
    
    % Trim the signal from the detected starting point
    signal = signal(startIndex:end);

    % Parse input arguments
    doPlot  = p.Results.doPlot;
    doTable = p.Results.doTable;
    
    % Parse vectors of input levels Lsk and Lnk
    Lsk = p.Results.Lsk(:).';
    Lnk = p.Results.Lnk(:).';
    
    adjustAmbientNoiseFlag = false;
    adjustAuditoryMaskingFlag = false;
    
    if any(~isnan(Lsk)) && any(~isnan(Lnk)) % Both Lsk and Lnk given
        adjustAmbientNoiseFlag = true;
        adjustAuditoryMaskingFlag = true;
        Isk = 10 .^ (Lsk / 10);
        Ink = 10 .^ (Lnk / 10);
    elseif any(~isnan(Lsk)) && any(isnan(Lnk)) % Only Lsk given
        adjustAuditoryMaskingFlag = true;
        Isk = 10 .^ (Lsk / 10);
        Ink = zeros(1, 7);
    elseif any(isnan(Lsk)) && any(~isnan(Lnk)) % Only Lnk given
        warning("Ambient noise levels alone (Lnk), without signal levels" + ...
            " (Lsk), are insufficient for calculating the ambient noise" + ...
            " adjustment. Therefore, the adjustment step will be omitted.")
    end

    octaveBands_all = [125, 250, 500, 1000, 2000, 4000, 8000]; % octave bands in Hz
    fm_all = [1.6, 1, 0.63, 2, 1.25, 0.8, 2.5; ... % modulation frequencies in Hz
            8, 5, 3.15, 10, 6.25, 4, 12.5];
    
    % Band-filter the input signal and cut the first 200 ms to suppress the
    % transient effects of the used IIR octave filters
    signalFiltered = bandFiltering(signal, fs, octaveBands_all);
    signalFiltered = signalFiltered(0.2 * fs : end, :);
    
    % Detect the envelope
    signalEnvelope = envelopeDetection(signalFiltered, fs);

    % Compute modulation depths of the input signal
    mk_o = MTF(signalEnvelope, fs, fm_all);
    
    % Compute modulation depths of the reference signal if it was passed to the STIPA function
    if ~isnan(p.Results.reference)
        reference = p.Results.reference;
        fsRef     = p.Results.fsRef;
        
        referenceFiltered = bandFiltering(reference, fsRef, octaveBands_all);
        referenceEnvelope = envelopeDetection(referenceFiltered, fsRef);
        mk_i = MTF(referenceEnvelope, fsRef, fm_all);
        mk = mk_o ./ mk_i;
    else % Use the default modulation depth 0.55
        mk = mk_o ./ 0.55;
    end
    
    % Check for any value in 'mk' exceeding the threshold of 1.3
    if any(mk(:) > 1.3)
        warning(['One or more m-values are higher than 1.3, which is very ' ...
            'unlikely and suggests non-sinusoidal fluctuations or impulsive' ...
            ' noises. This indicates an invalid measurement!'])
    end
    
    % Adjust mk values to ambient noise
    if adjustAmbientNoiseFlag == true
        mk = adjustAmbientNoise(mk, Isk, Ink);
    end
    
    % Adjust mk values to auditory masking and threshold effects
    if adjustAuditoryMaskingFlag == true
        mk = adjustAuditoryMasking(mk, Lsk, Isk, Ink);
    end
    
    % Limit the value of modulation transfer values mk to avoid complex values in SNR
    % (Note that, in contrast to the paper, the limitation is applied after the adjustment steps)
    mk(mk > 1) = 1;
    
    % Calculate SNR from the modulation transfer values and limit the range to [-15; 15] dB
    SNR = computeSNR(mk);
    SNR = clipSNR(SNR);
    
    % Calculate Trasmission Index
    TI = computeTI(SNR);
    
    % Calculate Modulation Transmission Index
    MTI = computeMTI(TI);
    
    % Calculate the final value of Speech Transmission index
    STI = computeSTI(MTI, alpha_k, beta_k);

    % Plot and table display
    if doTable == 1
        displayTableSTI(mk, octaveBands_all, fm_all, MTI, STI, alpha_k, beta_k);
    end
    
    if doPlot == 1

        % Optional calibration: prompt user to select reference signal
        [calFile, calPath] = uigetfile({'*.wav;*.flac;*.mat','Audio/Mat files (*.wav,*.flac,*.mat)'}, ...
                                       'Select calibration signal (94 or 114 dB SPL). Cancel to skip.');
        if isequal(calFile,0)
            useCalibration = false;
            fprintf('No calibration signal provided. Levels will be relative to input units.\n');
        else
            useCalibration = true;
            % Load calibration data
            [calData, fsCal] = audioread(fullfile(calPath, calFile));
            % Resample if necessary
            if fsCal ~= fs
                warning('Calibration samplerate differs (fsCal=%d, fs=%d). Resampling...', fsCal, fs);
                calData = resample(calData, fs, fsCal);
            end
            % Prompt for calibration level
            calibLevels = {'94','114'};
            [sel, ok] = listdlg('PromptString','Select calibration level (dB SPL):', ...
                                'SelectionMode','single','ListString',calibLevels);
            if ok
                calLevel = str2double(calibLevels{sel});
            else
                calLevel = 94; % default
                warning('No level selected, assuming 94 dB SPL.\n');
            end
            % Compute reference pressure
            p_rms = 20e-6 * 10^(calLevel/20);
            % Compute calibration constant
            K = p_rms / rms(calData);
            % Apply calibration to entire signal
            signal = K * signal;   % now in Pa
            fprintf('Calibration applied using %s at %d dB SPL.\n', calFile, calLevel);
        end
    
        % Band filtering and transient removal
        signalFiltered = bandFiltering(signal, fs, octaveBands_all);
        signalFiltered = signalFiltered(0.2*fs:end, :);
    
        % Initialize level arrays
        Level  = zeros(1, size(signalFiltered,2));    % relative or SPL
        LevelA = zeros(1, size(signalFiltered,2));    % A-weighted levels

        % A-weighting correction
        A_factors = [-16.1, -8.6, -3.2, 0, 1.2, 1.0, -1.1];
    
        % Compute octave-band levels
        for k = 1:size(signalFiltered,2)
            prms = rms(signalFiltered(:,k));
            if useCalibration
                % Absolute level in dB SPL
                Level(k) = 20 * log10(prms / 20e-6);
            else
                % Relative level in dB (relative to input unit)
                Level(k) = 20 * log10(prms);
            end
            % A-weighted
            LevelA(k) = Level(k) + A_factors(k);
        end

        % Calculation of additional values if masking is applied
        if adjustAuditoryMaskingFlag
            % Compute La using helper function (only for the first 6 bands)
            La = computeLa(Lsk);
            a = 10.^(La/10);
            Iam = [0, Isk(1:end-1) .* a];       % Masking from the lower band
            Ak = [46, 27, 12, 6.5, 7.5, 8, 12]; % Threshold values (dB SPL)
            Irt = 10.^(Ak/10);                  % Threshold level
            Isum = I + Iam + Irt;               % Sum (received signal + noise)
            
            % For A-weighted levels – for the first 6 bands, extend to 7 elements (7th element = NaN)
            MaskingA = 10 * log10(Iam(2:end)) + A_factors(1:6);
            MaskingA = [MaskingA, NaN];
            ThresholdA = 10 * log10(Irt) + A_factors;
            TotalA   = 10 * log10(Isum) + A_factors;
        else
            % If masking is not applied, define empty values
            Iam = [];
            Irt = [];
            Isum = [];
            MaskingA = [];
            ThresholdA = [];
            TotalA = [];
        end
        
            % Call the displayPlotSTI function to display the graphs.
            displayPlotSTI(mk, signal, fs, ...
                          STI, MTI, Level, LevelA, ...
                          adjustAuditoryMaskingFlag, Iam, Irt, Isum, ...
                          MaskingA, ThresholdA, TotalA);
    end
end
%% Helper Functions

    function startIndex = detectSignalStart(sig, fs, threshold)
    % DETECTSIGNALSTART Detects the start of the signal based on 
    % an energy threshold.
        windowSize = round(0.1 * fs); % 100 ms window
        energy = movmean(sig.^2, windowSize);
        idx = find(energy > threshold, 1, 'first');
        if isempty(idx)
            startIndex = 1;
        else
            startIndex = idx;
        end
    end

    function y = bandFiltering(x, fs, octaveBands)
    % BANDFILTERING Filter input signal using a octave filter of 18th order 
    % to achieve a minimum of 42 dB attenuation at the center frequency 
    % of each adjacent band.
        
        filterOrder = 18;

        y = NaN(length(x), length(octaveBands));
        
        for bandIdx = 1:length(octaveBands)
            octFilt = octaveFilter(octaveBands(bandIdx), '1 octave', ...
                'SampleRate', fs, 'FilterOrder', filterOrder);
            y(:, bandIdx) = octFilt(x);
        end
    end

    function envelope = envelopeDetection(x, fs)
    % [ page 49 - A.3.1.2 - Annex A ]
    % ENVELOPEDETECTION Compute the intensity envelope by squaring 
    % the outputs of the bandpass filters and applying a low pass filter
    % at a cut-off frequency of 100 Hz.
        
        envelope = x .* x;
        envelope = lowpass(envelope, 100, fs);
    end

    function mk = MTF(signalEnvelope, fs, fm)
    % [ page 50 - A.3.1.3 - Annex A ]
    % MTF Calculate the modulation depths of the received signal's envelope 
    % for each octave band and modulation frequency.
        
        seconds = length(signalEnvelope) / fs; % duration of the signal in seconds
        
        mk = NaN(2, 7);
        
        for k = 1:7 % iterate over octave bands
            Ik = signalEnvelope(:, k); % signal envelope of k-th octave band
            
            for n = 1:2 % iterate over each modulation frequency in k-th octave band
                % Calculate the duration and index of the signal for a whole number of periods
                secondsWholePeriod = floor(fm(n, k) * seconds) / fm(n, k);
                indexWholePeriod = round(secondsWholePeriod * fs);
                t = linspace(0, secondsWholePeriod, indexWholePeriod).';
                
                % Calculate the modulation depths using a whole number of
                % periods for each specific modulation frequency
                mk(n, k) = 2 * sqrt(sum(Ik(1:indexWholePeriod) .* sin(2 * pi * fm(n, k) * t)) ^ 2 + ...
                           sum(Ik(1:indexWholePeriod) .* cos(2 * pi * fm(n, k) * t)) ^ 2) / sum(Ik(1:indexWholePeriod));
            end
            
        end
        
    end

    function mk_ = adjustAmbientNoise(mk, Isk, Ink)
    % [ page 48 - A.2.3 - Annex A ]
    % ADJUSTAMBIENTNOISE Adjust the m-values for ambient noise.   

        mk_ = mk .* (Isk ./ (Isk + Ink));
    end

    function mk_ = adjustAuditoryMasking(mk, Lsk, Isk, Ink)
    % ADJUSTAUDITORYMASKING Adjust the m-values for auditory masking 
    % and threshold effects.

        % [ page 48 - A.2.4 - Annex A ]
        Ik = Isk + Ink; % total acoustic intensity
       
        % Auditory masking as a function of the octave band level
        % [ page 53 - Table A.2 - Annex A ]
        La = computeLa(Lsk);

        % [ page 53 - A.4.2 - Annex A ]
        a = 10 .^ (La / 10);

        % Total acoustic intensity for the level-dependent auditory masking effect
        % [ page 52 - A.4.2 - Annex A ]
        Iamk = [0, Isk(1:end-1) .* a];
        
        % Absolute speech reception threshold
        % [ page 54 - Table A.3 - Annex A ]
        Ak = [46, 27, 12, 6.5, 7.5, 8, 12];

        % Acoustic intensity level of the reception threshold
        % [ page 54 - A.4.3 - Annex A ]
        Irtk = 10 .^ (Ak / 10);
        
        % [ page 48 - A.2.4 - Annex A ]
        mk_ = mk .* (Ik ./ (Ik + Iamk + Irtk));
    end

    function SNR = computeSNR(mk)
    % [ page 47 - A.2.1 - Annex A ]
    % COMPUTESNR Compute the Signal-to-Noise Ratio (SNR) from the MTF matrix
    % consisting of modulation ratios mk.
        
        SNR = 10 * log10(mk ./ (1 - mk));
    end

    function SNR_clipped = clipSNR(SNR)
    % [ page 47 - A.2.1 - Annex A ]
    % CLIPSNR Limit the SNR values to fit the range from -15 to 15 dB.
        
        SNR_clipped = SNR;
        SNR_clipped(SNR_clipped > 15) = 15;
        SNR_clipped(SNR_clipped < -15) = -15;
    end

    function TI = computeTI(SNR)
    % [ page 47 - A.2.1 - Annex A ]    
    % COMPUTETI Compute the Transmission index from the SNR.
        
        TI = (SNR + 15) / 30;
    end

    function MTI = computeMTI(TI)
    % [ page 47 - A.2.1 - Annex A ]
    % COMPUTEMTI Compute the Modulation Transfer index (MTI) from 
    % the Transmission index TI.

        MTI = mean(TI);
    end

    function La_vec = computeLa(L_vec)
    % [ page 53 - Table A.2 - Annex A ]
    % COMPUTELA Compute the La values for the first 6 octave bands.
        La_vec = zeros(1,6);
        for i = 1:6
            if L_vec(i) < 63
                La_vec(i) = 0.5 * L_vec(i) - 65;
            elseif L_vec(i) < 67
                La_vec(i) = 1.8 * L_vec(i) - 146.9;
            elseif L_vec(i) < 100
                La_vec(i) = 0.5 * L_vec(i) - 59.8;
            else
                La_vec(i) = -10;
            end
        end
    end

    function STI = computeSTI(MTI, alpha_k, beta_k)
    % [ page 46 - A.2.1 - Annex A ]
    % COMPUTESTI Compute the final Speech Transmission Index (STI) from 
    % the Modulation transfer indices MTI.
        STI = min(sum(alpha_k .* MTI) - sum(beta_k .* sqrt(MTI(1:end - 1) .* MTI(2:end))), 1);
    end