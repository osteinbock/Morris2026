% Morris and Steinbock, 2026, Figure 2c components plus aaditional graphs
% Computes clean time maps from experimental data and
% analyzes age data t_age along pattern edge and for the entire pattern

clear; clc; close all;

%% Define file parameters for all three datasets
fileStructs(1) = struct(...
    'fname', 'c3pt0M.xlsx', ...
    'label', '3 mol/L', ...
    'crop', true, ...
    'xmin', 1097, 'xmax', 2346, ...
    'ymin', 201,  'ymax', 1450);

fileStructs(2) = struct(...
    'fname', 'c5pt25M.xlsx', ...
    'label', '5.25 mol/L', ...
    'crop', true, ...
    'xmin', 897, 'xmax', 1946, ...
    'ymin', 401,  'ymax', 1450);

minArea = 15000;    % Minimum connected region area (in pixels)
scale = 0.028794;   % mm per pixel (5904 pix = 17.0 cm)

%% Loop through all datasets
for i = 1:length(fileStructs)
    
    fileStruct = fileStructs(i);
    
    %% Read Excel file data
    T = readtable(fileStruct.fname);
    data = table2array(T);
    
    %% Crop data if requested
    if fileStruct.crop
        data = data(fileStruct.ymin:fileStruct.ymax, fileStruct.xmin:fileStruct.xmax);
    end
    
    %% Remove small blobs
    processedData = data;
    binaryMask = (processedData ~= -10);
    cc = bwconncomp(binaryMask);
    stats = regionprops(cc, 'Area');
    for j = 1:length(stats)
        if stats(j).Area < minArea
            processedData(cc.PixelIdxList{j}) = -10;
        end
    end
    
    %% Create outer edge map (4-connectivity for thicker edge)
    validMask = (processedData ~= -10);
    edgeMask = bwperim(validMask, 4);
    
    % Create edge data: keep only edge pixels, set rest to -10
    edgeData = -10 * ones(size(processedData));
    edgeData(edgeMask) = processedData(edgeMask);
    
    %% Apply time reversal transformation to heatmaps
    % For processedData
    validMask1 = (processedData ~= -10);
    tmax1 = max(processedData(validMask1));
    tmin1 = min(processedData(validMask1));
    processedDataFlipped = processedData;
    processedDataFlipped(validMask1) = tmax1 - processedData(validMask1) + tmin1;
    
    % For edgeData
    validMask2 = (edgeData ~= -10);
    tmax2 = max(edgeData(validMask2));
    tmin2 = min(edgeData(validMask2));
    edgeDataFlipped = edgeData;
    edgeDataFlipped(validMask2) = tmax2 - edgeData(validMask2) + tmin2;
    
    %% Create display data with NaN for background (for gray background)
    processedDisplayFlipped = processedDataFlipped;
    processedDisplayFlipped(processedDataFlipped == -10) = NaN;
    
    %% Create scaled axes
    xVec = (1:size(processedData,2)) * scale;  % x-axis in mm
    yVec = (1:size(processedData,1)) * scale;  % y-axis in mm
    
    %% Prepare distribution data
    % For cleaned heatmap (inset)
    C1 = processedData(:);
    C1_valid = C1(C1 > 0);
    C1_sorted = sort(C1_valid);
    C1_min = min(C1_sorted);
    C1_flipped = max(C1_sorted) - C1_sorted + C1_min;
    % C1_flipped = C1_flipped(C1_flipped > C1_min);
    C1_flipped = flipud(C1_flipped);  % sort ascending
    x1 = (1:length(C1_flipped))';
    p1 = polyfit(x1, C1_flipped, 1);
    
    % For edge data (main right plot)
    C2 = edgeData(:);
    C2_valid = C2(C2 > 0);
    C2_sorted = sort(C2_valid);
    C2_min = min(C2_sorted);
    C2_flipped = max(C2_sorted) - C2_sorted + C2_min;
    % C2_flipped = C2_flipped(C2_flipped > C2_min);
    C2_flipped = flipud(C2_flipped);  % sort ascending
    x2 = (1:length(C2_flipped))';
    p2 = polyfit(x2, C2_flipped, 1);
    
    %% Publication-quality figure
    figure('Color', 'w', 'Units', 'centimeters', 'Position', [2 2 18 8], ...
           'Name', fileStruct.label);
    
    % Right panel: Edge distribution (create first to get reference size)
    ax2 = axes('Position', [0.58 0.15 0.38 0.75]);
    plot(x2, C2_flipped, 'k.', 'MarkerSize', 8);
    hold on;
    plot(x2, polyval(p2, x2), 'r-', 'LineWidth', 2);
    hold off;
    xlabel('\bfindex {\itn}', 'FontSize', 11);
    hy2 = ylabel('\bf{\itt}_{edge} (s)', 'FontSize', 11);
    hy2.Position(1) = hy2.Position(1) + 50;  % move rightward
    axis square;
    ylim([0 100]);
    set(gca, 'FontSize', 10, 'LineWidth', 1, 'Box', 'on');
    drawnow;
    ax2Pos = ax2.Position;
    
    % Left panel: Cleaned heatmap with gray background
    ax1 = axes('Position', [0.02 ax2Pos(2) ax2Pos(3)+0.08 ax2Pos(4)]);
    h1 = imagesc(xVec, yVec, processedDisplayFlipped);
    set(h1, 'AlphaData', ~isnan(processedDisplayFlipped));
    set(gca, 'Color', [0.7 0.7 0.7]);
    axis image;
    colormap('jet');
    clim([0 95]);
    xlabel('\bf{\itx} (mm)', 'FontSize', 12);
    ylabel('\bf{\ity} (mm)', 'FontSize', 12);
    set(gca, 'FontSize', 10, 'LineWidth', 1, 'Box', 'on');
    drawnow;
    ax1Pos = ax1.Position;
    
    % Add colorbar after and position manually
    cb1 = colorbar;
    cb1.Label.String = '\bf{\itt} (s)';
    cb1.Label.Interpreter = 'tex';
    cb1.Label.FontSize = 11;  % add this line
    cb1.Label.Rotation = 270;
    cb1.Label.VerticalAlignment = 'bottom';
    cb1.Label.Position(1) = 2.5;  % move label closer to colorbar
    cb1.Position = [ax1Pos(1)+ax1Pos(3)-0.05, ax1Pos(2), 0.02, ax1Pos(4)];
    set(ax1, 'Position', ax1Pos);  % restore ax1 position
    
    % Inset: Full distribution (upper left corner of ax2)
    insetWidth = 0.29;
    insetHeight = 0.29;
    ax3 = axes('Position', [ax2Pos(1)-0.04, ax2Pos(2)+ax2Pos(4)-insetHeight-0.035, insetWidth, insetHeight]);
    plot(x1/1e4, C1_flipped, 'k.', 'MarkerSize', 8);
    hold on;
    plot(x1/1e4, polyval(p1, x1), 'r-', 'LineWidth', 1.0);
    hold off;
    set(gca, 'FontSize', 8, 'LineWidth', 0.75, 'Box', 'on');
    set(gca, 'YAxisLocation', 'right');  % move y-axis to right
    axis square;
    hx3 = xlabel('\bf{\itn} × 10^4', 'FontSize', 10);
    hy3 = ylabel('\bf{\itt}_{all} (s)', 'FontSize', 10);
    set(hx3, 'Interpreter', 'tex');
    set(hy3, 'Interpreter', 'tex');
    ylim([0 100]);
    drawnow;  % force rendering before next iteration
    
    % Add subplot labels (a, b) outside boxes
    annotation('textbox', [ax1Pos(1)-0.02, ax1Pos(2)+ax1Pos(4), 0.05, 0.05], ...
        'String', 'a', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'EdgeColor', 'none', 'VerticalAlignment', 'bottom');
    annotation('textbox', [ax2Pos(1)-0.06, ax2Pos(2)+ax2Pos(4), 0.05, 0.05], ...
        'String', 'b', 'FontSize', 15, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'EdgeColor', 'none', 'VerticalAlignment', 'bottom');
        
end