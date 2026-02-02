% CA-based recreation of experimental tw-t heatmap
% Sweep over beta = 0.25, 0.5, 0.75 with fixed k*beta = 0.6
% Figure 1: 3x12 pattern grid (rows = beta, cols = seeds)
% Figure 2: 3x2 (left: heatmaps, right: age distributions)

clear all
close all

% === PARAMETERS ===
KB_PRODUCT = 1.0;
BETA_vals = [0.5, 0.35, 0.2]; %[0.75, 0.5, 0.25];
nBeta = length(BETA_vals);

% CA parameters
N = 250;
DELTA = 0;
BMIN = 2;
ATTCH = 0;
r0 = 10;

T_pregrow = 1000;
tw_vals = [15, 30, 45, 60, 90, 105, 135, 180, 225];  % 3x original
N_samples = 500;

% Time scaling for display (1 timestep = 3 seconds)
timeScale = 3;

% Seeds for 12 runs
seeds = [123, 456, 789, 101, 201, 303, 877, 321, 654, 443, 7, 870];
nRuns = length(seeds);

% Bin settings (3x original)
binEdges = 0:30:300;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
M = numel(binCenters);
N_tw = numel(tw_vals);

% === STORAGE ===
allPatterns = cell(nBeta, nRuns);
allA = cell(nBeta, nRuns);
allB = cell(nBeta, nRuns);
binnedDataAll = zeros(M, N_tw, nRuns, nBeta);
CtargetAll = cell(nBeta, nRuns);

% === CREATE FIGURE 1 BEFORE MAIN LOOP ===
fig1 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 14 4]);

marginL = 0.06;
marginR = 0.01;
marginT = 0.12;
marginB = 0.02;
gapH = 0.005;
gapV = 0.02;

panelW_fig1 = (1 - marginL - marginR - (nRuns-1)*gapH) / nRuns;
panelH_fig1 = (1 - marginT - marginB - (nBeta-1)*gapV) / nBeta;

% === MAIN LOOP OVER BETA VALUES ===
for bi = 1:nBeta
    BETA = BETA_vals(bi);
    KAGE = KB_PRODUCT / BETA;
    
    fprintf('\n##############################\n');
    fprintf('BETA = %.2f, k = %.2f\n', BETA, KAGE);
    fprintf('##############################\n');
    
    for ri = 1:nRuns
        seed = seeds(ri);
        fprintf('\n=== RUN %d/%d (seed=%d) ===\n', ri, nRuns, seed);
        
        fprintf('Pre-growing pattern for %d timesteps...\n', T_pregrow);
        
        rng(seed);
        
        A = ones(N,N);  newA = A;
        B = zeros(N,N); newB = B;
        C = zeros(N,N); newC = C;
        
        for y = 1:N
            for x = 1:N
                [~, rho] = cart2pol(x-N/2, y-N/2);
                if rho <= r0
                    A(y,x) = 0;
                    B(y,x) = 1;
                end
            end
        end
        
        for t = 1:T_pregrow
            if mod(t, 100) == 0
                fprintf('  t = %d/%d\n', t, T_pregrow);
            end
            
            for y = 2:N-1
                for x = 2:N-1
                    newA(y,x) = A(y,x);
                    newB(y,x) = B(y,x);
                    newC(y,x) = C(y,x);
                    sumA = sum(sum(A(y-1:y+1, x-1:x+1)));
                    sumB = sum(sum(B(y-1:y+1, x-1:x+1)));
                    if sumA > 0 && sumB > 0 && C(y,x) == 0
                        newA(y,x) = 0;
                        newB(y,x) = 0;
                        newC(y,x) = 1;
                    end
                    if C(y,x) > 0
                        newA(y,x) = 0;
                        newB(y,x) = 0;
                        newC(y,x) = C(y,x) + 1;
                    end
                end
            end
            
            xtarget = []; ytarget = []; Ctarget = [];
            m = 0;
            for y = 2:N-1
                for x = 2:N-1
                    if C(y,x) > 0 && sum(sum(B(y-1:y+1, x-1:x+1))) >= BMIN
                        m = m + 1;
                        xtarget(m) = x;
                        ytarget(m) = y;
                        Ctarget(m) = C(y,x);
                    end
                end
            end
            
            if m > 0
                Cprobs = exp(-KAGE * (Ctarget-1.0).^BETA);
                probs = Cprobs / sum(Cprobs);
                xtmp = mnrnd(1, probs);
                [~, ti] = max(xtmp);
                newB(ytarget(ti), xtarget(ti)) = 1;
                newC(ytarget(ti), xtarget(ti)) = 0;
                
                dsqfinal = inf;
                for y = 2:N-1
                    for x = 2:N-1
                        if A(y,x) == 1
                            dsq = sqrt((x-xtarget(ti))^2 + (y-ytarget(ti))^2);
                            if dsq < dsqfinal
                                dsqfinal = dsq;
                            end
                        end
                    end
                end
                
                xp = []; yp = []; p = 0;
                for y = 2:N-1
                    for x = 2:N-1
                        if A(y,x) == 1 && sum(sum(C(y-1:y+1, x-1:x+1))) >= ATTCH
                            dsq = sqrt((x-xtarget(ti))^2 + (y-ytarget(ti))^2);
                            if dsq <= dsqfinal + DELTA
                                p = p + 1;
                                xp(p) = x;
                                yp(p) = y;
                            end
                        end
                    end
                end
                
                ppick = randi([1 p]);
                newA(yp(ppick), xp(ppick)) = 0;
                newB(yp(ppick), xp(ppick)) = 0;
                newC(yp(ppick), xp(ppick)) = 1;
            end
            
            A = newA; B = newB; C = newC;
        end
        
        allPatterns{bi, ri} = C;
        allA{bi, ri} = A;
        allB{bi, ri} = B;
        
        fprintf('Pre-growth complete.\n');
        
        xtarget = []; ytarget = []; Ctarget = [];
        m = 0;
        for y = 2:N-1
            for x = 2:N-1
                if C(y,x) > 0 && sum(sum(B(y-1:y+1, x-1:x+1))) >= BMIN
                    m = m + 1;
                    xtarget(m) = x;
                    ytarget(m) = y;
                    Ctarget(m) = C(y,x);
                end
            end
        end
        
        CtargetAll{bi, ri} = Ctarget;
        fprintf('Found %d breach candidates.\n', m);
        
        binnedData = zeros(M, N_tw);
        
        for ti = 1:N_tw
            tw = tw_vals(ti);
            agedCtarget = Ctarget + tw;
            Cprobs = exp(-KAGE * (agedCtarget - 1.0).^BETA);
            probs = Cprobs / sum(Cprobs);
            
            breachAges = zeros(1, N_samples);
            for s = 1:N_samples
                xtmp = mnrnd(1, probs);
                [~, idx] = max(xtmp);
                breachAges(s) = Ctarget(idx);
            end
            
            binnedData(:, ti) = histcounts(breachAges, binEdges)';
        end
        
        binnedDataAll(:,:,ri,bi) = binnedData;
        
        fprintf('Min candidate age: %d, Max: %d\n', min(Ctarget), max(Ctarget));
        
        figure(fig1);
        
        left = marginL + (ri-1)*(panelW_fig1 + gapH);
        bottom = 1 - marginT - bi*panelH_fig1 - (bi-1)*gapV;
        
        ax = axes('Position', [left, bottom, panelW_fig1, panelH_fig1]);
        
        minC = min(C(C>0));
        maxC = max(C(C>0));
        dum1 = zeros(size(C));
        dum1(C>0) = 1 - (C(C>0)-minC) / (maxC-minC);
        img1 = cat(3, 0.8*A+dum1, 0.8*A+dum1*0.6+0.6*B, 0.8*A+B);
        
        imagesc(img1);
        axis equal off;
        
        if ri == 1
            text(-0.08, 0.5, sprintf('\\beta = %.2f', BETA_vals(bi)), ...
                'Units', 'normalized', ...
                'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial', ...
                'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
        
        drawnow;
        fprintf('  [%d/%d] Figure 1 updated\n', (bi-1)*nRuns + ri, nBeta*nRuns);
    end
end

fprintf('\n=== ALL RUNS COMPLETE ===\n');

%% === FIGURE 1: Finalize and Save ===
fprintf('Finalizing Figure 1...\n');

figure(fig1);
sgtitle(sprintf('All patterns (k\\beta = %.1f, %d seeds)', KB_PRODUCT, nRuns), ...
    'FontSize', 12, 'FontWeight', 'bold');

% print(fig1, 'Fig_CA_patterns_3x12.png', '-dpng', '-r300');
% fprintf('Saved: Fig_CA_patterns_3x12.png\n');

%% === FIGURE 2: 3x2 Heatmaps and Age Distributions (Publication Quality) ===
fprintf('Creating Figure 2: Heatmaps and age distributions...\n');

fig2 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 5.6 7]);
tiledlayout(3, 2, 'TileSpacing', 'tight', 'Padding', 'compact');

fontName = 'Arial';
fontSize = 10;
labelSize = 11;
gamma = 0.8;

for bi = 1:nBeta
    BETA = BETA_vals(bi);
    
    % --- Left column: Cumulative heatmap ---
    nexttile;
    
    binnedCumulative = sum(binnedDataAll(:,:,:,bi), 3);
    pCumulative = binnedCumulative ./ sum(binnedCumulative, 1);
    
    imagesc(pCumulative.^gamma);
    colormap(gca, turbo);
    caxis([0 1]);
    cb = colorbar;
    cb.Label.String = '{\bf\it P}';
    cb.Label.Interpreter = 'tex';
    cb.FontSize = fontSize - 1;
    
    yTickLabel = arrayfun(@(i)sprintf('%d–%d', binEdges(i), binEdges(i+1)), 1:M, 'UniformOutput', false);
    
    set(gca, 'XTick', 1:2:N_tw, 'XTickLabel', tw_vals(1:2:end));
    set(gca, 'YTick', 1:2:M, 'YTickLabel', yTickLabel(1:2:end));
    
    axis('square')
    xlabel('{\bf{\it t_w} (time steps)}', 'Interpreter', 'tex', 'FontSize', labelSize);
    ylabel('{\bf{\it t} (time steps)}', 'Interpreter', 'tex', 'FontSize', labelSize);
    
    text(0.05, 0.08, sprintf('\\beta = %.2f', BETA), ...
        'Units', 'normalized', 'FontSize', fontSize, 'FontWeight', 'bold', ...
        'FontName', fontName, 'Color', 'w');
    
    panelLabelsLeft = {'a', 'b', 'c'};
    text(-0.5, 1.05, panelLabelsLeft{bi}, 'Units', 'normalized', ...
        'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    % --- Right column: Age distributions ---
    nexttile;
    hold on;
    
    allSorted = cell(1, nRuns);
    maxLen = 0;
    minLen = inf;
    maxC = 0;
    for rj = 1:nRuns
        Csorted = sort(CtargetAll{bi, rj}(:));
        allSorted{rj} = Csorted;
        maxLen = max(maxLen, length(Csorted));
        minLen = min(minLen, length(Csorted));
        maxC = max(maxC, max(Csorted));
    end
    
    allData = zeros(nRuns, minLen);
    for rj = 1:nRuns
        allData(rj,:) = allSorted{rj}(1:minLen);
    end
    meanVals = mean(allData, 1);
    stdVals = std(allData, 0, 1);
    
    fill([1:minLen, fliplr(1:minLen)], [meanVals-stdVals, fliplr(meanVals+stdVals)], ...
        [0.8 0.2 0.2], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    plot(1:minLen, meanVals, '-', 'Color', [0.7 0.1 0.1], 'LineWidth', 1.5);
    
    plot([0 maxLen], [0 maxLen], 'k--', 'LineWidth', 1.2);
    
    xlim([0 400]);
    ylim([0 800]);
    axis square;
    
    set(gca, 'FontName', fontName, 'FontSize', fontSize);
    
    box on
    xlabel('{\bf index \it{n}}', 'Interpreter', 'tex', 'FontSize', labelSize);
    ylabel('{\bf{sorted \it C}}', 'Interpreter', 'tex', 'FontSize', labelSize);  

    panelLabelsRight = {'d', 'e', 'f'};
    text(-0.38, 1.05, panelLabelsRight{bi}, 'Units', 'normalized', ...
        'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    hold off; 
    
    drawnow;
    fprintf('  Figure 2: beta = %.2f complete\n', BETA);
end

% print(fig2, 'Fig_CA_heatmaps_3x2_kB1p0LONG_b.png', '-dpng', '-r300');
% print(fig2, 'Fig_CA_heatmaps_3x2_kB1p0LONG_b.pdf', '-dpdf', '-vector');
% fprintf('Saved: Fig_CA_heatmaps_3x2_kB1p0LONG_b.png and .pdf\n');