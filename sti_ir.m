function [STI, mk, STI_PA] = sti_ir(IR, fs, varargin)
% STI_IR calculates the Speech Transmission Index (STI) from a system's 
% impulse response (IR). The calculation follows the indirect (IR-based)
% method as specified in the IEC 60268-16 (ed.3) standard.
%
%   [STI, MK] = STI_IR(IR, FS) filters the IR into standard octave bands
%   and modulation frequencies, computes the modulation transfer function
%   (MTF) matrix MK via Schroeder’s method, and returns the overall STI.
%
%   [STI, MK] = STI_IR(IR, FS, 'Lsk', LSK) computes the Speech
%   Transmission Index with adjustment of the MTF for auditory masking and 
%   threshold effects.
%
%   [STI, MK] = STI_IR(IR, FS, 'Lsk', LSK, 'Lnk', LNK) computes the
%   Speech Transmission Index with adjustments of the MTF for ambient 
%   noise, and auditory masking and threshold effects.
%
%   [STI, MK] = STI_IR(IR, FS, 'doPlot', DOPLOT) when 'doPlot' is set to 1,
%   the function displays graphical output of the results.
%   In this mode, the user is prompted to select an optional calibration
%   reference signal (94 or 114 dB SPL) to compute absolute dB SPL levels.
%   If the prompt is canceled, levels are shown relative to the input units.
%   When doPlot is set to 0 (default), the graphical output is suppressed.
%
%   [STI, MK] = STI_IR(IR, FS, 'doTable', DOTABLE) displays 
%   the MTF matrix in a table when DOTABLE is set to 1 . 
%   If DOTABLE is set to 0 (default), the table is not shown.
%
% The computation follows the IEC 60268‑16 specification for indirect
% STI calculation from room impulse responses.
%
% Reference:
%   EN IEC 60268‑16:2020 Sound system equipment — Part 16: Objective rating of speech
%   intelligibility by speech transmission index.
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025

    % Check number of input arguments
    narginchk(2, 10);

    % Parse input arguments:
    p = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
    validColumnVector = @(x) isnumeric(x) && iscolumn(x);
    valid7PosVector   = @(x) isnumeric(x) && length(x) == 7 && all(x > 0);
    validBinaryFlag   = @(x) isnumeric(x) && isscalar(x) && (x==0 || x==1);

    addRequired(p, 'IR', validColumnVector);
    addRequired(p, 'fs', validScalarPosNum);
    addParameter(p, 'Lsk', NaN, valid7PosVector);
    addParameter(p, 'Lnk', NaN, valid7PosVector);
    addParameter(p, 'doTable', 0, validBinaryFlag);
    addParameter(p, 'doPlot',  0, validBinaryFlag);
    parse(p, IR, fs, varargin{:});

    % Definition of weighting and redundancy factors
    % weighting factors for male speech according to the standard
    % [ page 48 - Annex A ]
    alpha_k = [0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173];
    
    % redundancy factors for male speech according to the standard
    % [ page 48 - Annex A ]
    beta_k = [0.085, 0.078, 0.065, 0.011, 0.047, 0.095];

    % Parse input arguments
    doTable = p.Results.doTable;
    doPlot  = p.Results.doPlot;

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

% Zero-pad the IR so that the MTF computation covers the entire IR duration
IR = [IR; zeros(1.6 * fs, 1)];
len = size(IR,1);

% Time in seconds for each sample
time=((1:len)-1)'./fs;

% List of modulation frequencies
fm = [0.63, 0.8, 1, 1.25, 1.6, 2, 2.5, 3.15, 4, 5, 6.3, 8, 10, 12.5];
% fm = 10.^((-2:11)/10); % exact frequencies

% Octave band nominal center frequencies
octaveBands = [125, 250, 500, 1000, 2000, 4000, 8000];

filterOrder = 18; % in-band filter pseudo-order

% Apply octave band filtering using the integrated subfunction
octFilters = octbandfilter(IR, fs, octaveBands, filterOrder);

% Number of whole number cycles to use for each modulation frequency
fm_cycles = floor(len .* fm./fs);
% Number of samples to use for each modulation frequency
fm_len = floor(fs.*fm_cycles ./ fm);

% % Initialize the resulting MTF matrix
%     mk = zeros(length(fm), length(octaveBands));
%     for j = 1:length(fm)
%         % Compute the MTF using Schroeder's formula:
%         % Numerator: Sum over time of the squared filtered signal weighted by
%         % an exponential phase correction.
%         num = abs( sum( octFilters(1:fm_len(j), 1, :) .^2 .* ...
%                         exp(-2i*pi*fm(j)*time(1:fm_len(j)) ) ) );
%         den = sum( octFilters(1:fm_len(j), 1, :) .^2 );
%         mk(j, :) = num ./ den;
%         % octFilters(...) .^2  => represents p^2(t) in the given octave band
%         % exp(-2i*pi*fm(j)*time(1:fm_len(j))) => represents e^(-j 2π f_m t)
%         % sum(...) => the discrete analog of the integral ∫ p^2(t) e^(-j2πf_m t) dt
%         % abs(...) => takes the magnitude (absolute value) of the complex sum
%         % mk(j, :) = num ./ den => normalizes by the total energy, ∫ p^2(t) dt
%     end

    % Initialize the resulting MTF matrix
    mk = zeros(length(fm), length(octaveBands));
    
    for j = 1:length(fm)
        for k = 1:length(octaveBands)
            % p^2(t) for this modulation frequency and octave band k
            p2 = octFilters(1:fm_len(j), 1, k).^2;
    
            % Summation for sin and cos
            num_sin = sum( p2 .* sin(2*pi*fm(j)*time(1:fm_len(j))) );
            num_cos = sum( p2 .* cos(2*pi*fm(j)*time(1:fm_len(j))) );
    
            % Numerator
            % NOTE: In the indirect (IR-based) STI calculation method according to 
            % IEC 60268-16, the factor of 2 in the numerator is not applied. This method 
            % uses the "Schroeder sum" of p^2(t) and implicitly assumes a normalized 
            % input modulation. The factor of 2 would be used if dealing with a truly 
            % 100% modulated test signal and an intensity (peak-to-peak) modulation 
            % definition, which applies to the direct STI method rather than the 
            % indirect IR-based method.

            num = sqrt( num_sin^2 + num_cos^2 );
    
            % Denominator
            den = sum( p2 );
    
            % Final MTF
            mk(j, k) = num / den;
        end
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

    idx_pairs = [
      5, 12;   % 1.6  Hz & 8    Hz → band 125 Hz
      3, 10;   % 1    Hz & 5    Hz → band 250 Hz
      1,  8;   % 0.63 Hz & 3.15 Hz → band 500 Hz
      6, 13;   % 2    Hz & 10   Hz → band 1 kHz
      4, 11;   % 1.25 Hz & 6.3  Hz → band 2 kHz
      2,  9;   % 0.8  Hz & 4    Hz → band 4 kHz
      7, 14    % 2.5  Hz & 12.5 Hz → band 8 kHz
    ];

    MTI_STIPA = zeros(1,7);
    for k = 1:7
      MTI_STIPA(k) = mean([ TI(idx_pairs(k,1), k), TI(idx_pairs(k,2), k) ]);
    end
    
    % Calculate the final value of Speech Transmission index
    % [ page 46 - A.2.1 - Annex A ]
    STI = computeSTI(MTI, alpha_k, beta_k);

    STI_PA = sum(alpha_k .* MTI_STIPA) ...
          - sum(beta_k .* sqrt(MTI_STIPA(1:end-1) .* MTI_STIPA(2:end)));

    % Plot and table display
    if doTable == 1
        displayTableSTI(mk, octaveBands, fm, MTI, STI, alpha_k, beta_k, STI_PA);
    end
    if doPlot == 1
        % --- Optional calibration: prompt user to select reference signal ---
        [calFile, calPath] = uigetfile({'*.wav;*.flac;*.mat','Audio/Mat files (*.wav,*.flac,*.mat)'}, ...
                                       'Select calibration file (94 or 114 dB SPL). Cancel to skip.');
        if isequal(calFile,0)
            useCalibration = false;
            fprintf('No calibration signal provided. Levels will be relative to input units.\n');
        else
            useCalibration = true;
            [calData, fsCal] = audioread(fullfile(calPath, calFile));
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
                calLevel = 94;
                warning('No level selected, assuming 94 dB SPL.\n');
            end
            % Compute reference RMS pressure and calibration constant
            p_rms = 20e-6 * 10^(calLevel/20);
            K = p_rms / rms(calData);
            % Apply calibration to IR
            IR = K * IR;   % now in Pa
            fprintf('Calibration applied using %s at %d dB SPL.\n', calFile, calLevel);
        end

        % 1) Filter (calibrated) IR into octave bands
        P_oct = octbandfilter(IR, fs, octaveBands, filterOrder);
        P_oct = P_oct(ceil(0.2*fs):end, 1, :);

        % 2) Compute band energies
        I = squeeze(mean(P_oct.^2, 1))';     % linear energies
        if useCalibration
            Level = 10*log10(I / (20e-6)^2); % convert to dB SPL
        else
            Level = 10*log10(I);             % relative dB to input units
        end
        A_factors = [-16.1, -8.6, -3.2, 0, 1.2, 1.0, -1.1];
        LevelA = Level + A_factors;          % SPL dB(A)

        % 3) Compute masking & noise terms if needed
        if adjustAuditoryMaskingFlag
            La   = computeLa(Lsk);
            a    = 10.^(La./10);
            % masking from the lower band
            Iam  = [0, Isk(1:end-1).*a];
            % absolute threshold values
            Ak   = [46,27,12,6.5,7.5,8,12];
            Irt  = 10.^(Ak./10);            % threshold PSD
            Isum = I + Iam + Irt;           % total S+N+masking

            % A-weighted versions
            MaskingA   = [10*log10(Iam(2:end)) + A_factors(1:6), NaN];
            ThresholdA = 10*log10(Irt)    + A_factors;
            TotalA     = 10*log10(Isum)   + A_factors;
        else
            % no masking → pass empty arrays
            Iam = []; Irt = []; Isum = [];
            MaskingA=[]; ThresholdA=[]; TotalA=[];
        end

        % 4) Call the displayPlotSTI function to display the graphs.
        displayPlotSTI(mk, IR, fs, STI, MTI, Level, LevelA, ...
                      adjustAuditoryMaskingFlag, Iam, Irt, Isum, MaskingA, ThresholdA, TotalA);
    end
end
%% Helper Functions

    function P_octave = octbandfilter(IR, fs, octaveBands, filterOrder)
    % OCTBANDFILTER Filters the input signal IR into octave bands 
    % using octaveFilter.
        n_samples = length(IR);
        n_bands = length(octaveBands);
        P_octave = zeros(n_samples, 1, n_bands);
        
        for k = 1:n_bands
            % Create an octave filter using the specified filter order
            octFilt = octaveFilter(octaveBands(k), '1/2 octave', 'SampleRate', fs, 'FilterOrder', filterOrder);
            P_octave(:, 1, k) = octFilt(IR);
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