% Morris and Steinbock, 2026, Figures 3a-c
% Includes breach-formation stats 'bft' for all six concentrations c
% Produces stretched-exponential fits and summarizes results as
% beta(c), k(c) plots with square root fits as well as
% k(beta) (inset) with hyperbola fit

clear; clc; close all;

% 1M
bft5_1M = [1,1,1,1,1,4,2,12,0,1,1,1,5,1,1,0,0,1,1,0,1,1,0,0];
bft10_1M = [1,0,0,2,7,2,3,6,2,2,2,0,0,1,0,1,0,0,1,1,0,1,1];
bft15_1M = [1,2,2,14,2,2,2,0,1,0,1,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,1,0];
bft20_1M =[1,1,0,1,1,0,0,0,0,120,18,0,54,0,1,2,1,1];
bft30_1M = [1,2,2,1,2,1,1,40,0,0,0,0,1,1,0,1,1,1,1,1,1,1];
bft35_1M = [1,1,1,0,1,0,0,1,0,1,1,1,1,0,42,0,0,0,19,0,2,1,1];
bft45_1M = [2,2,22,2,13,0,1,0,1,0,0,1,0,1,0,4,1,15,1,0,0,1,0,19];
bft60_1M = [1,2,1,0,13,0,0,0,7,0,18,0,3,2,0,1,0,29,1,1,1,20];
bft75_1M = [2,1,1,1,2,52,1,0,63,16,0,0,39,86,0,0,74,0,0,0,1,0,89,55,8,0,0,0,16,45];

% 3M
bft5_3M = [0,0,1,63,0,0,5,0,0,1,1,1,0,0,1,0,1,72,2,0,0,2];
bft10_3M = [0,1,1,0,0,5,5,0,16,0,0,1,1,1,0,1,1,1,0,44,0,1];
bft15_3M = [0,0,21,27,27,56,0,0,0,1,0,0,0,0,0,0,0,40,0,1,8,1];
bft20_3M =[0,0,0,3,1,11,23,7,0,9,0,0,9,0,0,1,4,1,1,1,1,1];
bft30_3M = [0,0,7,0,0,9,24,1,17,0,2,1,1,1,1,1];
bft35_3M = [2,86,0,1,0,6,1,2,1,0,1,0,0,1,0,1,1,1,2,1,2,2];
bft45_3M = [1,32,5,7,0,0,10,0,0,30,8,1,0,0,1,0,1,0,1,0,6,75];
bft60_3M = [1,0,9,1,0,8,56,11,1,65,1,0,14,1,1,1,48,0,1,1,1,55,1];
bft75_3M = [43,90,77,46,39,46,21,0,0,1,30,33,0,0,65];

% 4M 
bft5_4M = [0,0,0,6,0,3,1,0,4,1,1,12,8,1,1,7,1,1,1,1,1];
bft10_4M = [0,8,0,1,0,0,0,0,1,1,0,10,1,3,5,3,1,0,15,1,3];
bft15_4M = [0,0,0,0,1,0,0,1,13,1,7,1,7,1,0,3,6,16,1,1,1,0,5];
bft20_4M = [0,1,65,1,1,1,1,1,0,0,8,1,4,3,1,8,49,4,1,1,1];
bft30_4M = [1,1,1,8,1,1,11,1,86,20,15,1,0,1,0,1,3];
bft35_4M =[0,1,22,1,10,4,1,1,1,1,1,77,1,0,1,9,5,25,9,0];
bft45_4M = [7,8,7,0,1,4,72,7,27,13,4,1,1,0,74,0,1];
bft60_4M = [62,1,42,1,0,1,1,34,15,1,62,55,0,12,1,6,1,0,21];
bft75_4M = [33,27,1,21,33,4,0,2,20,9,47,68,2,1,1,1,77,4,0,69,1];

% 4.5M
bft5_4p5M = [8,0,1,27,1,0,10,3,1,1,73,0,1,1,1,1,0,1,2,13,1];
bft10_4p5M = [0,10,3,1,0,22,67,1,1,1,1,1,6,1,1];
bft15_4p5M = [1,7,4,0,76,8,6,16,1,1,1,27,1,83,1,0,1,16,1,4,7,3,34,6];
bft20_4p5M =[0,1,14,9,12,1,1,1,4,5,5,74,22,0,6,1,29];
bft30_4p5M = [51,63,1,0,2,0,0,0,9,68,1,1,21,12,0,3,13,1,0,60,8];
bft35_4p5M = [12,12,1,1,0,1,46,17,0,3,1,7,4,0,1,26,4,73,4,1];
bft45_4p5M = [16,17,0,1,1,1,1,1,1,4,1,67,1,1,1,52,1,14];
bft60_4p5M = [31,45,4,15,1,1,74,0,0,1,1,81,0,17,1,3,33,0];
bft75_4p5M = [19,41,76,0,71,79,79,15,83,38,73,16,85,47,1,44,1,0,22,49,80,32];

% 5M
bft5_5M = [1,1,1,1,1,1,19,20,2,2,70,5,1,2,1,3,27,1,4,6,4,12,1,82,9,1,6,11,49,18,58,48,8,21,1,8,7,2,1,9,3,19,65];
bft10_5M = [34,80,35,2,5,4,33,15,2,2,16,20,1,79,2,3,36,5,28,65,2,2,1,1,14,92,1,9,1,3,14,20,2,7,1,48,31,1,54,28,15,57,1,3,4,1,39,2,55,1,3,1,8,7,2,1,9,1,1,21,75,40];
bft15_5M = [21,1,12,3,46,48,33,39,57,7,3,13,20,12,22,18,57,29,30,20,4,16,16,52,13,52,1,32,9,54,54,45,17,45,1,13,2,10,20,14,60,51,29,33,46,46,17,21,27,1,34,3,16,11,1,5,90,3,34,17,74,1,54,1,2,62,58,74,36,67,3,1];
bft20_5M = [73,4,19,15,55,24,50,15,1,25,39,50,58,13,1,1,47,40,28,48,31,1,6,15,11,16,30,29,19,1,8,1,1,43,55,54,22,2,10,29,2,13,6,14,41,1,33,32,12,2,12,4,4,27];
bft30_5M = [11,1,4,56,56,47,3,3,15,13,47,49,29,11,75,24,55,1,5,54,1,16,23,78,75,14,48,23,2,8,2,1,1,73,51,64,42,1,1,14,42,78,32,38,33,4,21,33,32,50,1,1,50];
bft35_5M =[60,47,67,13,35,9,45,3,17,1,18,9,3,44,5,7,10,11,46,4,58,12,51,63,76,46,4,11,0,0,25,4,9,51,28,12,51,2,13,5,91,73,16,1,74,25,17,2,21,2];
bft45_5M = [3,12,3,1,27,7,36,2,50,80,1,1,35,94,73,39,3,37,1,65,18,17,11,1,25,6,14,16,41,62,1,42,30,71,28,8,20,1,1,13,44,52,9,87,66,14,47,1,72,89,1,7,14,45,53,54,29,1];
bft60_5M = [1,1,2,23,42,74,29,75,45,46,83,13,3,18,20,4,39,46,11,1,5,66,95,7,18,24,86,41,20,1,44,93,62,19,42,2,36,1,16,1,51];
bft75_5M = [5,5,66,10,10,19,3,22,77,71,86,91,64,11,27,61,3,72,2,41,58,1,29,31,90,8,1,31,45,77,46,42,80,21,3,4,14,7,67,45,3,5,39,4,19,4,54];

% 5.25M
bft5_5p25M = [1,1,1,1,0,3,1,0,1,1,0,1,0,1,7,0,1,1,0,2,1,1,1,1,5];
bft10_5p25M = [0,1,1,1,3,18,0,1,80,32,1,1,0,1,3,1,6,0];
bft15_5p25M = [1,1,21,0,1,0,59,3,1,13,1,16,1,0,74,0,8,3];
bft20_5p25M =[1,7,15,1,38,0,17,0,0,15,33,1,1,23,37];
bft30_5p25M = [0,23,19,35,1,1,55,1,0,16,7,1,0,75,0,0,0,20];
bft35_5p25M = [0,0,0,47,75,1,1,1,1,13,9,45,1,0,1,7,0,1,2,5,1,40,50];
bft45_5p25M = [10,1,1,1,84,0,38,38,17,81,1,12,1,0,1,26,1,93,0,1,14,27,2,1,93,10,87,3,10,5];
bft60_5p25M = [57,63,76,45,1,1,4,47,0,22,44,1,70,1,28,37,3,34,54,25,1,31,0,1,1,0,9];
bft75_5p25M = [7,75,1,70,1,78,12,84,91,56,1,26,47,17,0,0,47,45,11,26,61,60,15,66,44,1,4,83,33,42,65,49,68];

timeVals = [5,10,15,20,30,35,45,60,75];

%% ---- Organize data by concentration ----
allBFT_1M = {bft5_1M, bft10_1M, bft15_1M, bft20_1M, bft30_1M, bft35_1M, bft45_1M, bft60_1M, bft75_1M};
allBFT_3M = {bft5_3M, bft10_3M, bft15_3M, bft20_3M, bft30_3M, bft35_3M, bft45_3M, bft60_3M, bft75_3M};
allBFT_4M = {bft5_4M, bft10_4M, bft15_4M, bft20_4M, bft30_4M, bft35_4M, bft45_4M, bft60_4M, bft75_4M};
allBFT_4p5M = {bft5_4p5M, bft10_4p5M, bft15_4p5M, bft20_4p5M, bft30_4p5M, bft35_4p5M, bft45_4p5M, bft60_4p5M, bft75_4p5M};
allBFT_5M = {bft5_5M, bft10_5M, bft15_5M, bft20_5M, bft30_5M, bft35_5M, bft45_5M, bft60_5M, bft75_5M};
allBFT_5p25M = {bft5_5p25M, bft10_5p25M, bft15_5p25M, bft20_5p25M, bft30_5p25M, bft35_5p25M, bft45_5p25M, bft60_5p25M, bft75_5p25M};

concNames = {'1M', '3M', '4M', '4.5M', '5M', '5.25M'};
allData = {allBFT_1M, allBFT_3M, allBFT_4M, allBFT_4p5M, allBFT_5M, allBFT_5p25M};
cVals = [1, 3, 4, 4.5, 5, 5.25];

%% ---- bin settings ----
binEdges   = 0:10:100;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
M = numel(binCenters);
N = numel(timeVals);

%% ---- tick labels ----
xTick = 1:N;
xTickLabel = arrayfun(@num2str, timeVals, 'UniformOutput', false);
yTickLabel = arrayfun(@(i)sprintf('%d–%d',binEdges(i),binEdges(i+1)),1:M,'UniformOutput',false);

%% ---- fitting setup ----
t  = binCenters(:);
tw = timeVals(:)';
normCols = @(X) bsxfun(@rdivide, X, sum(X,1)+eps);
opts = statset('MaxIter',20000,'Display','off');

%% ---- Storage arrays ----
pExpAll = cell(1, numel(allData));
pStrExpFitAll = cell(1, numel(allData));
kStrAll = zeros(1, numel(allData));
betaStrAll = zeros(1, numel(allData));
rmseStrExpAll = zeros(1, numel(allData));

%% ---- Main fitting loop ----
for c = 1:numel(allData)
    allBFT = allData{c};
    
    % Bin & normalize
    binnedData = zeros(M, N);
    for j = 1:N
        binnedData(:,j) = histcounts(allBFT{j}, binEdges);
    end
    pExp = binnedData ./ sum(binnedData,1,'omitnan');
    pExpAll{c} = pExp;
    yData = pExp(:);
    
    % Stretched-exponential fit
    modStrExp = @(params) normCols(exp(-params(1)*(t+tw).^params(2)));
    modFunStrExp = @(params,~) reshape(modStrExp(params), [],1);
    initialGuess = [10.0, 0.5];
    paramsStr = nlinfit([],yData,modFunStrExp,initialGuess,opts);
    kStrAll(c) = paramsStr(1);
    betaStrAll(c) = paramsStr(2);
    pStrExpFitAll{c} = modStrExp(paramsStr);
    
    % Compute RMSE
    rmseStrExpAll(c) = sqrt(mean((yData - pStrExpFitAll{c}(:)).^2));
end

%% ---- Profile likelihood for beta (to get proper error bars) ----
fprintf('Computing profile likelihood for β...\n');
betaScan = linspace(0.02, 0.55, 100);
rmseScanAll = zeros(numel(allData), numel(betaScan));

for c = 1:numel(allData)
    for ib = 1:numel(betaScan)
        betaFixed = betaScan(ib);
        modK = @(k) normCols(exp(-k*(t+tw).^betaFixed));
        modFunK = @(k,~) reshape(modK(k), [],1);
        kOpt = nlinfit([], pExpAll{c}(:), modFunK, 10, opts);
        fitted = modK(kOpt);
        rmseScanAll(c,ib) = sqrt(mean((pExpAll{c}(:) - fitted(:)).^2));
    end
end

% Extract profile-likelihood based confidence intervals
betaProfileMin = zeros(1, numel(allData));
betaErrLo = zeros(1, numel(allData));
betaErrHi = zeros(1, numel(allData));
rmseThreshold = 1.05;  % 5% above minimum

for c = 1:numel(allData)
    [minRMSE, idx] = min(rmseScanAll(c,:));
    betaProfileMin(c) = betaScan(idx);
    
    % Find confidence interval (where RMSE < threshold * min)
    inCI = rmseScanAll(c,:) < rmseThreshold * minRMSE;
    ciIndices = find(inCI);
    if ~isempty(ciIndices)
        betaLo = betaScan(ciIndices(1));
        betaHi = betaScan(ciIndices(end));
        betaErrLo(c) = betaProfileMin(c) - betaLo;
        betaErrHi(c) = betaHi - betaProfileMin(c);
    end
end

% %% ---- Figure: All Heatmaps (Experimental vs Fit) ----
% figure('Color','w','Units','normalized','Position',[.05 .1 .9 .6]);
% 
% for c = 1:numel(allData)
%     % Experimental data
%     subplot(2,6,c);
%     imagesc(pExpAll{c}); axis tight square;
%     colormap(turbo); caxis([0 1]);
%     set(gca, ...
%         'XTick', xTick, 'XTickLabel', xTickLabel, ...
%         'YTick', 1:M, 'YTickLabel', yTickLabel, ...
%         'FontSize', 7);
%     xlabel('{\it t}_{\rm w} (s)','Interpreter','tex','FontSize',8);
%     ylabel('{\it t} (s)','Interpreter','tex','FontSize',8);
%     title(sprintf('%s - Data', concNames{c}), 'FontWeight', 'bold', 'FontSize', 9);
%     colorbar('FontSize',6);
% 
%     % Stretched-exponential fit
%     subplot(2,6,6+c);
%     imagesc(pStrExpFitAll{c}); axis tight square;
%     colormap(turbo); caxis([0 1]);
%     set(gca, ...
%         'XTick', xTick, 'XTickLabel', xTickLabel, ...
%         'YTick', 1:M, 'YTickLabel', yTickLabel, ...
%         'FontSize', 7);
%     xlabel('{\it t}_{\rm w} (s)','Interpreter','tex','FontSize',8);
%     ylabel('{\it t} (s)','Interpreter','tex','FontSize',8);
%     title(sprintf('%s - Fit (\\beta=%.2f)', concNames{c}, betaStrAll(c)), 'FontWeight', 'bold', 'FontSize', 9);
%     colorbar('FontSize',6);
% end
% 
% sgtitle('Breach Formation Time: P(t | t_w) \propto exp(-k(t+t_w)^\beta)', 'FontWeight', 'bold', 'FontSize', 14);
% 
% %% ---- Figure 2: Profile likelihood ----
% figure('Color','w','Units','normalized','Position',[.05 .1 .9 .35]);
% 
% for c = 1:numel(allData)
%     subplot(1,6,c);
%     plot(betaScan, rmseScanAll(c,:), 'b-', 'LineWidth', 1.5);
%     hold on;
%     [minRMSE, ~] = min(rmseScanAll(c,:));
%     yline(rmseThreshold * minRMSE, 'r--', 'LineWidth', 1);
%     xline(betaProfileMin(c), 'k--', 'LineWidth', 1);
%     xlabel('\beta');
%     ylabel('RMSE');
%     title(sprintf('%s', concNames{c}), 'FontWeight', 'bold');
%     axis square;
%     xlim([0 0.55]);
% end
% sgtitle('Profile Likelihood: RMSE vs \beta (k optimized at each \beta)', 'FontWeight', 'bold', 'FontSize', 12);

%% ---- Print summary to console ----
fprintf('\n========================================\n');
fprintf('STRETCHED-EXPONENTIAL FIT RESULTS\n');
fprintf('Model: P(t|t_w) ∝ exp(-k(t+t_w)^β)\n');
fprintf('========================================\n\n');

fprintf('β values (profile likelihood):\n');
fprintf('%-6s | %8s | %18s | %s\n', 'Conc', 'β', '95% CI', 'Note');
fprintf('-------|----------|--------------------|-----------\n');
for c = 1:numel(allData)
    note = '';
    if c == 5
        note = 'unreliable';
    end
    fprintf('%-6s | %8.3f | [%6.3f, %6.3f] | %s\n', concNames{c}, betaProfileMin(c), ...
        betaProfileMin(c)-betaErrLo(c), betaProfileMin(c)+betaErrHi(c), note);
end


%% ---- Figure: heatmaps (1M, 4.5M, 5M) ----
pubConcs = [1, 4, 5];  % indices for 1M, 4.5M, 5M
pubNames = {'1 M', '4.5 M', '5 M'};
figure('Color','w','Units','inches','Position',[1 1 8.5 4.5]);

% Adjust spacing
gapH = 0.0;   % horizontal gap between columns
gapV = 0.06;  % vertical gap between rows
marginL = 0.08;
marginR = 0.12;
marginT = 0.08;
marginB = 0.14;
panelW = (1 - marginL - marginR - 3*gapH) / 3;
panelH = (1 - marginT - marginB - gapV) / 2;

% Gamma for nonlinear color scaling (<1 expands low values)
gamma = 0.8;

for i = 1:3
    c = pubConcs(i);
    
    % Top row: Experimental data
    ax1 = axes('Position', [marginL + (i-1)*(panelW+gapH), marginB + panelH + gapV, panelW, panelH]);
    imagesc(pExpAll{c}.^gamma);
    axis tight equal;
    colormap(turbo); caxis([0 1]);
    set(gca, ...
        'XTick', 1:2:N, 'XTickLabel', xTickLabel(1:2:end), ...
        'YTick', 1:2:M, 'YTickLabel', yTickLabel(1:2:end), ...
        'FontSize', 10, 'FontName', 'Arial', ...
        'TickLength', [0.02 0.02], 'LineWidth', 0.75, 'Box', 'on');
    if i == 1
        yl = ylabel('{\it t} (s)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
        yl.Position(1) = yl.Position(1) - 0.3;
    else
        set(gca, 'YTickLabel', []);
    end
    set(gca, 'XTickLabel', []);
    title(pubNames{i}, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    % Bottom row: Stretched-exponential fit
    ax2 = axes('Position', [marginL + (i-1)*(panelW+gapH), marginB, panelW, panelH]);
    imagesc(pStrExpFitAll{c}.^gamma);
    axis tight equal;
    colormap(turbo); caxis([0 1]);
    set(gca, ...
        'XTick', 1:2:N, 'XTickLabel', xTickLabel(1:2:end), ...
        'YTick', 1:2:M, 'YTickLabel', yTickLabel(1:2:end), ...
        'FontSize', 10, 'FontName', 'Arial', ...
        'TickLength', [0.02 0.02], 'LineWidth', 0.75, 'Box', 'on');
    xlabel('{\it t}_w (s)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
    if i == 1
        yl = ylabel('{\it t} (s)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
        yl.Position(1) = yl.Position(1) - 0.3;
    else
        set(gca, 'YTickLabel', []);
    end
    
    % Add beta label
    text(0.95, 0.08, sprintf('\\beta = %.2f', betaStrAll(c)), ...
        'Units', 'normalized', 'FontSize', 11, 'FontName', 'Arial', ...
        'HorizontalAlignment', 'right', 'Color', 'w', 'FontWeight', 'bold');
end

% Single colorbar
cbAx = axes('Position', [1-marginR+0.02, marginB, 0.02, 2*panelH+gapV]);
Pvals = linspace(0, 1, 256)';
imagesc(cbAx, [0 1], [0 1], Pvals.^gamma);
set(cbAx, 'YDir', 'normal');
colormap(cbAx, turbo);
caxis(cbAx, [0 1]);
set(cbAx, 'XTick', [], ...
    'YTick', [0, 0.2, 0.4, 0.6, 0.8, 1], ...
    'YTickLabel', {'0', '0.2', '0.4', '0.6', '0.8', '1'}, ...
    'YAxisLocation', 'right', ...
    'FontSize', 10, 'FontName', 'Arial', ...
    'Box', 'on', ...
    'LineWidth', 0.75, ...
    'TickLength', [0 0]);  % removes all tick marks
ylabel(cbAx, '\it{P}', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');

% Export
% print(gcf, 'Fig4_heatmaps.png', '-dpng', '-r300');
% % print(gcf, 'Fig4_heatmaps.pdf', '-dpdf', '-vector');
% fprintf('Saved: Fig4_heatmaps.png (not Fig4_heatmaps.pdf)\n');


%% ---- Figure: β vs concentration ----
figure('Color','w','Units','inches','Position',[1 1 4.5 4]);

% Plot data with profile likelihood error bars
errorbar(cVals, betaProfileMin, betaErrLo, betaErrHi, 'ko', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 6, 'CapSize', 10);
hold on;

% Fit: beta = c1 * sqrt(c2 - c) excluding 5M
fitIdx = [1,2,3,4,6];  % exclude 5M
fitfun = @(p, c) p(1) * sqrt(p(2) - c);
p0 = [0.1, 6];
pFit = nlinfit(cVals(fitIdx), betaProfileMin(fitIdx), fitfun, p0);
cFit = linspace(0.5, pFit(2)-0.0001, 250);
betaFit = fitfun(pFit, cFit);
plot(cFit, betaFit, 'r-', 'LineWidth', 2);

% Formatting
set(gca, ...
    'FontSize', 12, 'FontName', 'Arial', ...
    'LineWidth', 1, 'Box', 'on', ...
    'TickLength', [0.02 0.02]);
xlabel('{\it c} (M)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
% ylabel('\beta', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');  % larger to appear bolder
ylabel('\fontsize{16}\bf\beta', 'FontName', 'Arial');

xlim([0 6]);
ylim([0 0.5]);
axis square;

% Add fit equation with proper square root
text(0.4, 0.05, '$\beta = 0.19\sqrt{5.6 - c/(1\,\mathrm{M})}$', ...
    'FontSize', 12, 'Interpreter', 'latex');
text(-0.185, 1.002, 'b', 'Units', 'normalized', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial');

% Export
% print(gcf, 'Fig5_beta_vs_c.png', '-dpng', '-r300');
% print(gcf, 'Fig5_beta_vs_c.pdf', '-dpdf', '-vector');
% fprintf('Saved: Fig5_beta_vs_c.png and Fig5_beta_vs_c.pdf\n');

% -----------------------------------------------------------------------

%% ---- Profile likelihood for k ----
fprintf('Computing profile likelihood for k...\n');
kScan = linspace(1, 40, 100);
rmseScanK = zeros(numel(allData), numel(kScan));

for c = 1:numel(allData)
    for ik = 1:numel(kScan)
        kFixed = kScan(ik);
        modBeta = @(beta) normCols(exp(-kFixed*(t+tw).^beta));
        modFunBeta = @(beta,~) reshape(modBeta(beta), [],1);
        betaOpt = nlinfit([], pExpAll{c}(:), modFunBeta, 0.3, opts);
        fitted = modBeta(betaOpt);
        rmseScanK(c,ik) = sqrt(mean((pExpAll{c}(:) - fitted(:)).^2));
    end
end

% Extract profile-likelihood based confidence intervals for k
kProfileMin = zeros(1, numel(allData));
kErrLo = zeros(1, numel(allData));
kErrHi = zeros(1, numel(allData));

for c = 1:numel(allData)
    [minRMSE, idx] = min(rmseScanK(c,:));
    kProfileMin(c) = kScan(idx);
    
    inCI = rmseScanK(c,:) < rmseThreshold * minRMSE;
    ciIndices = find(inCI);
    if ~isempty(ciIndices)
        kLo = kScan(ciIndices(1));
        kHi = kScan(ciIndices(end));
        kErrLo(c) = kProfileMin(c) - kLo;
        kErrHi(c) = kHi - kProfileMin(c);
    end
end

%% ---- Figure 6: Publication-quality k vs concentration with inset ----
figure('Color','w','Units','inches','Position',[1 1 4.5 4]);
% Main plot: k vs c with error bars
errorbar(cVals, kProfileMin, kErrLo, kErrHi, 'ko', ...
'MarkerFaceColor', 'k', 'LineWidth', 1.5, 'MarkerSize', 6, 'CapSize', 10);
hold on;
% Fit: k = k0 / sqrt(c* - c) with fixed exponent -1/2
fitIdx = [1,2,3,4,6];  % exclude 5M
cStar = pFit(2);
fitfun_k = @(k0, c) k0 ./ sqrt(cStar - c);
p0_k = 5;
k0_fit = nlinfit(cVals(fitIdx), kProfileMin(fitIdx), fitfun_k, p0_k);
cFit = linspace(0.5, cStar-0.3, 250);
kFit = fitfun_k(k0_fit, cFit);
plot(cFit, kFit, 'r-', 'LineWidth', 2);
% Formatting
set(gca, ...
'FontSize', 12, 'FontName', 'Arial', ...
'LineWidth', 1, 'Box', 'on', ...
'TickLength', [0.02 0.02]);
xlabel('{\it c} (M)', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
ylabel('{\it k}', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
xlim([0 6]);
ylim([0 25]);
axis square;
% Equation - lower right
text(0.5, 2.55, sprintf('$k = %.1f\\,/\\sqrt{%.1f - c/(1\\,\\mathrm{M})}$', k0_fit, cStar), ...
    'FontSize', 12, 'Interpreter', 'latex');
% Inset: β vs k (square, y-axis on right)
% axInset = axes('Position', [0.2, 0.57, 0.3, 0.3]);
axInset = axes('Position', [0.18, 0.565, 0.32, 0.32]);
plot(betaProfileMin, kProfileMin, 'ko', 'MarkerFaceColor', 'k', 'LineWidth', 1.2, 'MarkerSize', 5);
hold on;
% Label select points
text(betaProfileMin(1)-0.03, kProfileMin(1)-1.9, '1M', 'FontSize', 8);
text(betaProfileMin(5)-0.035, kProfileMin(5)+2.2, '5M', 'FontSize', 8);
set(axInset, ...
'FontSize', 9, 'FontName', 'Arial', ...
'LineWidth', 0.75, 'Box', 'on', ...
'TickLength', [0.02 0.02], ...
'YAxisLocation', 'right', ...
'XTick', [0, 0.2, 0.4], ...
'YTick', [0, 10, 20]);
xlabel('\beta', 'FontSize', 10, 'FontName', 'Arial');
ylabel('{\it k}', 'FontSize', 10, 'FontName', 'Arial');
xlim([0 0.5]);
ylim([0 25]);
axis square;
% In the inset section, after plotting the data points:
betaLine = linspace(0.05, 0.5, 100);
kbConst = mean(kProfileMin([1,2,3,4,6]) .* betaProfileMin([1,2,3,4,6]));  % ~2.8
kLine = kbConst ./ betaLine;
plot(betaLine, kLine, 'r-', 'LineWidth', 1);
text(0.24, 18, '$k\beta = 2.8$', 'FontSize', 10, 'Interpreter', 'latex');
text(-0.59, 1.12, 'c', 'Units', 'normalized', 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial');
% Export
% print(gcf, 'Fig6_k_vs_c.png', '-dpng', '-r300');
% print(gcf, 'Fig6_k_vs_c.pdf', '-dpdf', '-vector');
% fprintf('Saved: Fig6_k_vs_c.png and Fig6_k_vs_c.pdf\n');