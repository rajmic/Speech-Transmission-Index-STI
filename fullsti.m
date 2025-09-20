function [STI, mk] = fullsti(signal, fs, varargin)
% FULLSTI calculates the Full Speech Transmission Index (STI) based on
% a recorded signal. The calculation follows the IEC 60268-16 (ed.3) standard.
%
%   [STI, MK] = FULLSTI(SIGNAL, FS) computes the Speech Transmission Index
%   from the test signal SIGNAL and its sampling frequency FS in Hz,
%   returning the Speech Transmission index STI and a 14-by-7 matrix of the
%   respective modulation transfer values MK for each octave band and
%   modulation frequency.
%
%   [STI, MK] = FULLSTI(SIGNAL, FS, 'SegmentDuration', SEGMENTDURATION, ...
%   'SilenceDuration', SILENCEDURATION) sets the duration of each segment
%   and the silence between segments (in seconds). Default values are
%   10 seconds (segment) and 0 seconds (silence).
%
%   [STI, MK] = FULLSTI(SIGNAL, FS, REFERENCE) computes the Speech
%   Transmission Index using the REFERENCE signal instead of the default
%   value 1, which results from the FULL STI test signal. If the sampling
%   frequency of the REFERENCE is not specified, it is assumed that the
%   sampling frequency is the same as the sampling frequency of the test
%   SIGNAL FS.
%
%   [STI, MK] = FULLSTI(SIGNAL, FS, REFERENCE, FSREF) computes the Speech
%   Transmission Index using the REFERENCE signal and its sampling
%   frequency FSREF.
%   
%   [STI, MK] = FULLSTI(SIGNAL, FS, 'Lsk', LSK) computes the Speech
%   Transmission Index with adjustment of the MTF for auditory masking and 
%   threshold effects.
%
%   [STI, MK] = FULLSTI(SIGNAL, FS, 'Lsk', LSK, 'Lnk', LNK) computes the
%   Speech Transmission Index with adjustments of the MTF for ambient 
%   noise, and auditory masking and threshold effects.
%
%   [STI, MK] = FULLSTI(SIGNAL, FS, 'SignalStart', SIGNALSTART) manually
%   specifies the start of the test signal in the recording (in seconds).
%   If not provided, the start is automatically detected using an energy
%   threshold.
%
%   [STI, MK] = FULLSTI(SIGNAL, FS, 'ManualTrimTime', TRIMTIME) specifies how
%   many seconds to trim from both the beginning and end of each segment
%   before further processing. Default value is 0.5 seconds.
%
%   [STI, MK] = FULLSTI(SIGNAL, FS, 'PaddingDuration', PADDINGDURATION) specifies 
%   how many seconds of silence to append at the end of the signal to ensure 
%   sufficient data for the last segment. Default value is 1 second.
%
%   [STI, MK] = FULLSTI(SIGNAL, FS, 'doPlot', DOPLOT) when 'doPlot' is set to 1
%   (default), the function displays graphical output of the results.
%   In this mode, the user is prompted to select an optional calibration
%   reference signal (94 or 114 dB SPL) to compute absolute dB SPL levels.
%   If the prompt is canceled, levels are shown relative to the input units.
%   When doPlot is set to 0, the graphical output is suppressed.
%
%   [STI, MK] = FULLSTI(SIGNAL, FS, 'doTable', DOTABLE) displays 
%   the MTF matrix in a table when DOTABLE is set to 1 (default). 
%   If DOTABLE is set to 0, the table is not shown.
%
% Reference:
%   EN IEC 60268-16:2020 Sound system equipment — Part 16: Objective rating of speech
%   intelligibility by speech transmission index.
%
% The code is based on the Pavel Záviška implementation.
% Záviška, P., Rajmic, P., Schimmel, J. Matlab Implementation of STIPA
% (Speech Transmission Index for Public Address Systems). AES Europe, 2024.
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025

    % Check number of input arguments
    narginchk(2, 16);
    
    % Parse input arguments:
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    validColumnVector = @(x) isnumeric(x) && iscolumn(x);
    valid7PosVector   = @(x) isnumeric(x) && length(x) == 7 && all(x > 0);
    validNonNegativeScalar = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
    validBinaryFlag   = @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1);
    
    addRequired(p,  'signal', validColumnVector);
    addRequired(p,  'fs', validScalarPosNum);
    addOptional(p,  'reference', NaN, validColumnVector);
    addOptional(p,  'fsRef', fs, validScalarPosNum);
    addParameter(p, 'Lsk', NaN, valid7PosVector);
    addParameter(p, 'Lnk', NaN, valid7PosVector);
    addParameter(p, 'SegmentDuration', 10, validScalarPosNum); 
    addParameter(p, 'SilenceDuration', 0, validScalarPosNum);
    addParameter(p, 'doPlot', 1, validBinaryFlag);
    addParameter(p, 'doTable', 1, validBinaryFlag);
    addParameter(p, 'SignalStart', NaN, validNonNegativeScalar);
    addParameter(p, 'ManualTrimTime', 0.5, validNonNegativeScalar);
    addParameter(p, 'PaddingDuration', 1, validScalarPosNum);
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
    
    % Append a padding of silence at the end based on user-defined PaddingDuration
    paddingDuration = p.Results.PaddingDuration; % in seconds (default = 1)
    paddingSamples = paddingDuration * fs;
    signal = [signal; zeros(paddingSamples, 1)];

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
    
    % Define segment and silence duration
    segmentDuration = p.Results.SegmentDuration; % user-provided or default segment duration in seconds
    silenceDuration = p.Results.SilenceDuration;  % user-provided or default silence duration in seconds
    segmentLength = segmentDuration * fs;  % segment length in sampled points
    silenceLength = silenceDuration * fs;  % silence length in sampled points
    
    octaveBands_all = [125, 250, 500, 1000, 2000, 4000, 8000]; % octave bands
    fm_all = [0.63, 0.8, 1, 1.25, 1.6, 2, 2.5, 3.15, 4, 5, 6.3, 8, 10, 12.5]; % modulation frequencies
    
    % Compute MTF matrix for the test signal using a dedicated function
    mk_o = computeMTFMatrix(signal, fs, segmentLength, silenceLength, octaveBands_all, fm_all, p.Results.ManualTrimTime);
    
    % Compute modulation depths for the reference signal if provided
    if ~isnan(p.Results.reference)
        reference = p.Results.reference;
        fsRef = p.Results.fsRef;
        segmentLengthRef = segmentDuration * fsRef;
        silenceLengthRef = silenceDuration * fsRef;
        mk_i = computeMTFMatrix(reference, fsRef, segmentLengthRef, silenceLengthRef, octaveBands_all, fm_all, p.Results.ManualTrimTime);
        mk = mk_o ./ mk_i;
    else
        mk = mk_o ./ 1;
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
    % [ page 47 - A.2.1 - Annex A ]
    SNR = computeSNR(mk);
    SNR = clipSNR(SNR);
    
    % Calculate Trasmission Index
    % [ page 47 - A.2.1 - Annex A ]
    TI = computeTI(SNR);
    
    % Calculate Modulation Transmission Index
    % [ page 47 - A.2.1 - Annex A ]
    MTI = computeMTI(TI);
    
    % Calculate the final value of Speech Transmission index
    % [ page 46 - A.2.1 - Annex A ]
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
        signalFiltered = bandFiltering(signal, fs);
        signalFiltered = signalFiltered(0.2*fs:end, :);

        % A-weighting correction
        A_factors = [-16.1, -8.6, -3.2, 0, 1.2, 1.0, -1.1];

        % Compute octave-band levels
        Level = zeros(1,7);
        for k = 1:7
            prms = rms(signalFiltered(:,k));
            if useCalibration
                % Absolute level in dB SPL
                Level(k) = 20 * log10(prms / 20e-6);
            else
                % Relative level in dB (relative to input unit)
                Level(k) = 10 * log10(mean(signalFiltered(:,k).^2));
            end
        end
        % A-weighted
        LevelA = Level + A_factors;

        if adjustAuditoryMaskingFlag
            % Compute La using helper function (only for the first 6 bands)
            La = computeLa(Lsk);
            a = 10.^(La/10);
            Isk = 10.^(Lsk/10);
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
                      adjustAuditoryMaskingFlag, Iam, Irt, Isum, MaskingA, ThresholdA, TotalA);
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

    function mkMatrix = computeMTFMatrix(sig, fs, segmentLen, silenceLen, octaveBand, fm, manualTrimTime)
    % COMPUTEMTFMATRIX Computes the MTF matrix for a given signal.
    numFm = length(fm);
    numOct = length(octaveBand);

    % Initialize MTF matrix
    mkMatrix = NaN(numFm, numOct);

    filterOrder = 18;
    % Pre-generation of filters for all 7 octave bands
    octFilters = cell(1, numOct);
    for i = 1:numOct
         octFilters{i} = octaveFilter(octaveBand(i), '1/2 octave', 'SampleRate', fs, 'FilterOrder', filterOrder);
    end

    % Sequence for process each segment
    segIdx = 1; % Initial index for segments
        for fmIdx = 1:numFm
            for octIdx = 1:numOct
                % Determine the start and end of the segment
                segmentStart = (segIdx - 1) * (segmentLen + silenceLen) + 1;
                segmentEnd = segmentStart + segmentLen - 1;

                if segmentEnd > length(sig)
                warning(['Insufficient data for a complete segment. Stopping ' ...
                    'segmentation. Consider increasing the "PaddingDuration" ' ...
                    'parameter (e.g., set it to a value greater than 1) to ' ...
                    'ensure sufficient data for the final segment.']);
                    break;
                end

                % Extract the signal segment (each segment contains 'segmentDuration' seconds of signal)
                currentSegment = sig(segmentStart:segmentEnd);

                % Additional trimming: remove manualTrimTime seconds from the beginning and end of the segment
                if manualTrimTime > 0
                    trimSamples = round(manualTrimTime * fs);
                    if length(currentSegment) > 2 * trimSamples
                        currentSegment = currentSegment(trimSamples+1:end-trimSamples);
                    else
                        warning('Segment is too short for additional trimming. Skipping trimming.');
                    end
                end
                                               
                % Apply band filtering for the current segment
                % [ page 49 - A.3.1.2 - Annex A ]
                % Filter input signal using an octave filter for the given octave band
                % with a filter order of 18th to achieve a minimum of 42 dB attenuation 
                % at the center frequency of each adjacent band.
                signalFiltered = octFilters{octIdx}(currentSegment);
                
                % Remove first 200 ms to suppress filter transients
                signalFiltered = signalFiltered(0.2 * fs : end, :);
                
                % Envelope detection
                % [ page 49 - A.3.1.2 - Annex A ]
                signalEnvelope = envelopeDetection(signalFiltered, fs);
                
                % Compute MTF for the current segment
                % [ page 50 - A.3.1.3 - Annex A ]
                mkMatrix(fmIdx, octIdx) = MTF(signalEnvelope, fs, fm(fmIdx), octaveBand(octIdx));
                
                segIdx = segIdx + 1; % Move to the next segment
            end
        end
        % Warning when segment is shorter than 10 s
        trimmedDuration = length(currentSegment) / fs;
        if trimmedDuration < 10
            warning(['Segment duration after trimming is only %.2f s, which is less than ' ...
                     'the recommended 10 s (IEC 60268-16).'], trimmedDuration);
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

    function mk = MTF(signalEnvelope, fs, fm, octaveBands)
    % [ page 50 - A.3.1.3 - Annex A ]    
    % MTF Calculate the modulation depths of the received signal's envelope
    % for each octave band and modulation frequency.
                      
        seconds = length(signalEnvelope) / fs; % duration of the signal in seconds
        mk = NaN(length(fm), length(octaveBands)); % Initialize matrix to store modulation depths
       
        for k = 1:length(octaveBands) % iterate over octave bands
            Ik = signalEnvelope(:, k); % signal envelope of k-th octave band
            
            for n = 1:length(fm) % iterate over each modulation frequency in k-th octave band
                % Calculate the duration and index of the signal for a whole number of periods
                secondsWholePeriod = floor(fm(n) * seconds) / fm(n);
                indexWholePeriod = round(secondsWholePeriod * fs);
                t = linspace(0, secondsWholePeriod, indexWholePeriod).';
                
                % Calculate the modulation depths using a whole number of
                % periods for each specific modulation frequency
                modDepth = 2 * sqrt(sum(Ik(1:indexWholePeriod) .* sin(2 * pi * fm(n) * t)) .^ 2 + ...
                sum(Ik(1:indexWholePeriod) .* cos(2 * pi * fm(n) * t)) .^ 2) / sum(Ik(1:indexWholePeriod));
                
                mk(n,k) = modDepth; % Store the value in MTF matrix
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
        Irtk = 10 .^ (Ak ./ 10);
        
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

    function y = bandFiltering(x, fs)
    % BANDFILTERING Filter input signal using a octave filter of 18th order 
    % to achieve a minimum of 42 dB attenuation at the center frequency 
    % of each adjacent band.
        
        filterOrder = 18;
        octaveBands = [125, 250, 500, 1000, 2000, 4000, 8000];
        
        y = NaN(length(x), length(octaveBands));
        
        for bandIdx = 1:length(octaveBands)
            octFilt = octaveFilter(octaveBands(bandIdx), '1 octave', ...
                'SampleRate', fs, 'FilterOrder', filterOrder);
            y(:, bandIdx) = octFilt(x);
        end

    end

    function STI = computeSTI(MTI, alpha_k, beta_k)
    % [ page 46 - A.2.1 - Annex A ]
    % COMPUTESTI Compute the final Speech Transmission Index (STI) from 
    % the Modulation transfer indices MTI.
        STI = min(sum(alpha_k .* MTI) - sum(beta_k .* sqrt(MTI(1:end - 1) .* MTI(2:end))), 1);
    end