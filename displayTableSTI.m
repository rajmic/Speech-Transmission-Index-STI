function displayTableSTI(mkMatrix, octBands, fm, MTI, STI, alpha_k, beta_k, varargin)
% DISPLAYTABLE Displays MTF, MTI, weighting factors, STI, CIS, ALCons and category.
%   Supports Full STI, STIPA, and indirect IR-based methods (STI_IR).
%
% Usage:
%   displayTableSTI(mkMatrix, octBands, fm, MTI, STI, alpha_k, beta_k)
%   displayTableSTI(..., STIPA_IR)
%
% Inputs:
%   mkMatrix  – modulation transfer matrix
%                * FULL STI: size 14×7, rows correspond to fm vector (14×1)
%                * STIPA   : size 2×7,  rows correspond to low/high fm per band
%                * IR-based: size 14×7, rows correspond to fm vector (14×1) from IR method
%   octBands  – 1×7 vector of octave band center frequencies [125 250 … 8000]
%   fm        – modulation frequencies:
%                * FULL STI: 14×1 vector [0.63 0.8 … 12.5]
%                * STIPA   : 2×7 matrix, low/high modulation freqs per band
%   MTI       – 1×7 modulation transfer index per octave band
%   STI       – scalar speech transmission index
%   alpha_k   – 1×7 weighting factors for STI calculation
%   beta_k    – 1×6 (or padded 1×7) redundancy factors for STI
%   STIPA_IR  – (optional) scalar STIPA index from IR-based method
%
% Example:
%   % Full STI:
%   displayTableSTI(mk_full, [125 250 500 1000 2000 4000 8000], fm_full, MTI_full, STI_full, alpha_k, beta_k);
%   % STIPA:
%   displayTableSTI(mk_stipa, [125 250 500 1000 2000 4000 8000], fm_stipa, MTI_stipa, STI_stipa, alpha_k, beta_k);
%   % IR-based STI:
%   displayTableSTI(mk_ir,   [125 250 500 1000 2000 4000 8000], fm_ir,   MTI_ir,   STI_ir,   alpha_k, beta_k, STIPA_IR);
%
% Copyright Šimon Cieslar, Brno University of Technology, 2024-2025

    % Check optional STIPA_IR input
    hasSTIPA_IR = ~isempty(varargin);
    if hasSTIPA_IR
        STIPA_IR = varargin{1};
    end

    % Compute CIS (Common Intelligibility Scale) from STI
    CIS = 1 + log10(STI);

    % Compute Articulation Loss of Consonants (ALCons) in %
    ALCons = 170.5405 * exp(-5.419 * STI);

    % Determine STI category for normal listeners
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

    % Determine CIS category
    if CIS < 0.48
        cisCategory = 'bad';
    elseif CIS < 0.65
        cisCategory = 'poor';
    elseif CIS < 0.78
        cisCategory = 'fair';
    elseif CIS < 0.88
        cisCategory = 'good';
    else
        cisCategory = 'excellent';
    end

    % Determine table dimensions
    [nRows, nCols] = size(mkMatrix);

    % Column names from octave bands
    colNames = arrayfun(@(x) sprintf('%d Hz', x), octBands, 'UniformOutput', false);

    % Determine data and row names based on fm input and mkMatrix size
    if isvector(fm) && numel(fm) == nRows
        % FULL STI or IR-based: fm is vector of length nRows
        tableData1 = mkMatrix;
        rowNames1 = arrayfun(@(x) sprintf('%.2f Hz', x), fm, 'UniformOutput', false);
    elseif ismatrix(fm) && all(size(fm) == [nRows, nCols]) && nRows == 2
        % STIPA case: fm is 2×nCols matrix
        % Combine MTF and modulation frequencies into one table
        tableData1 = [mkMatrix; fm];
        rowNames1 = {'MTF values for each first mod freq', 'MTF values for each second mod freq', 'First modulation freq (Hz)', 'Second modulation freq (Hz)'};
    else
        % Fallback generic names
        tableData1 = mkMatrix;
        rowNames1 = arrayfun(@(i) sprintf('Row %d', i), 1:nRows, 'UniformOutput', false);
    end

    % Create figure
    fig = figure('Name','MTF & STI Details','Units','normalized','Position',[0.1 0.1 0.8 0.8]);

    % Title for MTF table
    uicontrol('Parent',fig,'Style','text','String','MTF per Modulation Frequency',...        
        'Units','normalized','Position',[0.01 0.95 0.98 0.04],'FontSize',12,'FontWeight','bold');

    % First table: MTF values (and fm for STIPA)
    uitable('Parent',fig,'Data',tableData1,'ColumnName',colNames,'RowName',rowNames1,...
        'Units','normalized','Position',[0.01 0.45 0.98 0.50],'ColumnWidth','auto');

    % Title for MTI and factors table
    uicontrol('Parent',fig,'Style','text','String','MTI & Weighting/Redundancy Factors',...        
        'Units','normalized','Position',[0.01 0.42 0.98 0.04],'FontSize',12,'FontWeight','bold');
                
    % Second table: MTI, alpha_k, beta_k
    beta_padded = [beta_k, NaN];
    colNames2 = arrayfun(@(x) [num2str(x) ' Hz'], octBands, 'UniformOutput', false);
    rowNames2 = {'MTI','alpha_k','beta_k'};
    tableData2 = [MTI; alpha_k; beta_padded];
    uitable('Parent',fig,'Data',tableData2,'ColumnName',colNames2,'RowName',rowNames2,...
        'Units','normalized','Position',[0.01 0.25 0.98 0.15],'ColumnWidth','auto');
        
    % STI, CIS, ALCons & Categories table
    uicontrol('Parent', fig, 'Style', 'text', 'String', 'STI, CIS, ALCons & Categories', ...
        'Units', 'normalized', 'Position', [0.01 0.20 0.98 0.04], 'FontSize', 12, 'FontWeight', 'bold');
    if hasSTIPA_IR
        colNames3 = {'STI', 'STIPA(IR)', 'CIS', 'ALCons (%)', 'STI Category', 'CIS Category'};
        data3 = {STI, STIPA_IR, CIS, ALCons, stiCategory, cisCategory};
    else
        colNames3 = {'STI', 'CIS', 'ALCons (%)', 'STI Category', 'CIS Category'};
        data3 = {STI, CIS, ALCons, stiCategory, cisCategory};
    end
    uitable('Parent', fig, 'Data', data3, 'ColumnName', colNames3, 'RowName', {' '}, ...
        'Units', 'normalized', 'Position', [0.01 0.05 0.98 0.10], 'ColumnWidth', 'auto');
end