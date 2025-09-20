function displayPlotSTI(mkMatrix, signal, fs, STI, MTI, Level, LevelA, adjustAuditoryMaskingFlag, Iam, Irt, Isum, MaskingA, ThresholdA, TotalA)
% DISPLAYPLOTSTI Displays waveform, MTF pixel map, STI, category, MTI bar chart, and SPL/A-weighted SPL plots.
%   Supports Full STI, STIPA and IR-based (indirect) STI layouts, with optional auditory masking details.
%
% Usage:
%   displayPlotSTI(mkMatrix, signal, fs, STI, MTI, Level, LevelA, false, [], [], [], [], [], [])
%   displayPlotSTI(mkMatrix, signal, fs, STI, MTI, Level, LevelA, true, Iam, Irt, Isum, MaskingA, ThresholdA, TotalA)
%   % For IR-based method, pass the impulse response as `signal` and include masking arrays:
%   displayPlotSTI(mkMatrix, IR, fs, STI, MTI, Level, LevelA, true, Iam, Irt, Isum, MaskingA, ThresholdA, TotalA)
%
% Inputs:
%   mkMatrix                  – modulation transfer matrix
%                                * FULL STI: size 14×7, rows correspond to fm vector (14×1)
%                                * STIPA   : size 2×7, rows correspond to low/high fm per band
%                                * IR-based: size 14×7, rows correspond to fm vector (14×1) from IR method
%   signal                    – input time-domain signal vector
%   fs                        – sampling frequency of  the signal (Hz)
%   STI                       – scalar speech transmission index
%   MTI                       – 1×7 modulation transfer index per octave band
%   Level                     – 1×7 sound pressure levels per band (dB)
%   LevelA                    – 1×7 A-weighted SPL per band (dB(A))
%   adjustAuditoryMaskingFlag – boolean flag to include auditory masking plots
%   Iam                       – 1×7 masking power spectral densities (linear)
%   Irt                       – 1×7 threshold power spectral densities (linear)
%   Isum                      – 1×7 total (signal+noise) power spectral densities (linear)
%   MaskingA                  – 1×7 A-weighted masking levels (dB(A))
%   ThresholdA                – 1×7 A-weighted threshold levels (dB(A))
%   TotalA                    – 1×7 A-weighted total S+N levels (dB(A))
%
% Example:
%   % Without masking details:
%   displayPlotSTI(mk_full, audioSignal, fs, STI_full, MTI_full, Level, LevelA, false);
%   % With masking details:
%   displayPlotSTI(mk_full, audioSignal, fs, STI_full, MTI_full, Level, LevelA, true, Iam, Irt, Isum, MaskingA, ThresholdA, TotalA);
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025

    % create figure in pixels
    fig = figure('Units','pixels');
    
    % get screen size: [x0 y0 totalWidth totalHeight]
    screenSize = get(0, 'ScreenSize');
    
    % desired window size
    winWidth  = 1200;
    winHeight = 800;
    
    % calculate bottom-left corner to center the window
    left   = (screenSize(3) - winWidth)  / 2;
    bottom = (screenSize(4) - winHeight) / 2;
    
    % set Position = [left bottom width height]
    fig.Position = [left, bottom, winWidth, winHeight];
    clf;

    % 1) Signal waveform
    subplot(3,2,1);
    plot((0:length(signal)-1)/fs, signal);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Signal waveform');
    grid on;

    % 2) Pixel map of MTF Values (dynamic orientation and annotations)
    figure(8)
    % subplot(3,2,2);
    [nRows, nCols] = size(mkMatrix);
    if nRows > nCols
        % FULL STI: transpose for wide layout
        dataImg = mkMatrix';
        xTicks = 1:nRows;
        yTicks = 1:nCols;
        xLabels = {'0.63','0.8','1','1.25','1.6','2','2.5','3.15','4','5','6.3','8','10','12.5'};
        yLabels = {'125','250','500','1k','2k','4k','8k'};
        xlabelTxt = 'Modulation Frequencies (Hz)';
        ylabelTxt = 'Octave Bands (Hz)';
        imagesc(dataImg);
        set(gca, 'XTick', xTicks, 'XTickLabel', xLabels);
        set(gca, 'YTick', yTicks, 'YTickLabel', yLabels);
    else
        % STIPA: direct layout with annotations inside cells
        dataImg = mkMatrix;
        imagesc(dataImg);
        % keep X axis for octave bands
        set(gca, 'XTick', 1:nCols, 'XTickLabel', {'125','250','500','1k','2k','4k','8k'});
        set(gca, 'YTick', []);
        xlabelTxt = 'Octave Bands (Hz)';
        ylabelTxt = 'Modulation Frequencies (Hz)';
        % annotate each cell with its modulation frequency label
        fmLabels = {
            '1.6','1','0.63','2','1.25','0.8','2.5';
            '8','5','3.15','10','6.25','4','12.5'};
        for r = 1:nRows
            for c = 1:nCols
                text(c, r, fmLabels{r,c}, 'HorizontalAlignment','center', ...
                    'VerticalAlignment','middle', 'FontSize',8, 'Color','k');
            end
        end
    end
    colormap(parula);
    colorbar;
    axis tight;
    title('MTF Values Pixel Map');
    xlabel(xlabelTxt);
    ylabel(ylabelTxt);
    clim([0 1]);

    figure(fig)

    % 3) Display STI value and category
    subplot(3,2,3);
    axis off;
    % Determine STI category label
    if STI < 0.30
        stiCategory = 'bad';
    elseif STI < 0.45
        stiCategory = 'poor';
    elseif STI < 0.60
        stiCategory = 'fair';
    elseif STI < 0.75
        stiCategory = 'good';
    else
        stiCategory = 'excellent';
    end
    % Display STI and category
    text(0.5, 0.55, sprintf('Speech Transmission Index (STI) = %.2f', STI), ...
         'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(0.5, 0.35, sprintf('Category: %s', stiCategory), ...
         'FontSize', 12, 'FontAngle', 'italic', 'HorizontalAlignment', 'center');

    % 4) Bar chart of MTI per octave band
    subplot(3,2,4);
    hBar = bar(MTI,'r');
    ylim([0 1]);
    set(gca,'XTick',1:7,'XTickLabel',{'125','250','500','1k','2k','4k','8k'});
    xlabel('Octave band (Hz)');
    ylabel('Modulation Transfer Index (MTI)');
    % Add value labels on top of each bar
    xt = hBar.XEndPoints;
    yt = hBar.YEndPoints;
    labels = arrayfun(@(v) sprintf('%.2f', v), MTI, 'UniformOutput', false);
    text(xt, yt, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',8);

    % 5) SPL plot
    subplot(3,2,5);
    hold on;
    grid on;
    plot(1:7, Level,'b-o','DisplayName','Received');
    if adjustAuditoryMaskingFlag
        plot(1:7,10*log10(Iam),'-x','DisplayName','Masking');
        plot(1:7,10*log10(Irt),':','DisplayName','Threshold');
        plot(1:7,10*log10(Isum),'--+','DisplayName','Total S+N');
    end
    hold off;
    set(gca,'XTick',1:7,'XTickLabel',{'125','250','500','1k','2k','4k','8k'});
    xlabel('Octave band (Hz)');
    ylabel('SPL (dB)');
    legend('Location','best');
    title('Sound Pressure Level (SPL)');

    % 6) A-weighted SPL plot
    subplot(3,2,6);
    hold on;
    grid on;
    plot(1:7, LevelA,'b-o','DisplayName','Received dB(A)');
    if adjustAuditoryMaskingFlag
        plot(1:7, MaskingA,'-x','DisplayName','Masking dB(A)');
        plot(1:7, ThresholdA,':','DisplayName','Threshold dB(A)');
        plot(1:7, TotalA,'--+','DisplayName','Total S+N dB(A)');
    end
    hold off;
    set(gca,'XTick',1:7,'XTickLabel',{'125','250','500','1k','2k','4k','8k'});
    xlabel('Octave band (Hz)');
    ylabel('SPL (dB(A))');
    legend('Location','best');
    title('A-weighted Sound Pressure Level (SPL)');
end