% Morris and Steinbock, 2026, Figure 4a
% CA: sweep different beta, k*beta pairs

clear all
close all

% === PARAMETER SETTINGS ===
KB_vals = [1.2, 0.6, 0.3];  % rows (top to bottom)
BETA_vals = [0.5, 0.25, 0.05];  % columns

nRows = length(KB_vals);
nCols = length(BETA_vals);

% === FIXED PARAMETERS ===
N = 150;
DELTA = 0;
BMIN = 2;
ATTCH = 0;
r0 = 10;
tMax = 40;

% === STORAGE ===
allPatterns = cell(nRows, nCols);
allA = cell(nRows, nCols);
allB = cell(nRows, nCols);

% === RUN ALL SIMULATIONS ===
for ri = 1:nRows
    KB_PRODUCT = KB_vals(ri);
    for ci = 1:nCols
        BETA = BETA_vals(ci);
        KAGE = KB_PRODUCT / BETA;
        
        fprintf('Row %d, Col %d: k*beta=%.2f, beta=%.2f, k=%.2f\n', ri, ci, KB_PRODUCT, BETA, KAGE);
        
        rng(123);
        
        A = ones(N,N);  newA = A;
        B = zeros(N,N); newB = B;
        C = zeros(N,N); newC = C;
        
        [X,Y] = meshgrid(1:N, 1:N);
        rC = (X-N/2).^2 + (Y-N/2).^2;
        rmaxC = 0;
        
        for y = 1:N
            for x = 1:N
                [~, rho] = cart2pol(x-N/2, y-N/2);
                if rho <= r0
                    A(y,x) = 0;
                    B(y,x) = 1;
                end
            end
        end
        
        t = 0;
        while rmaxC < (N/2-2)^2 % && t < tMax
            t = t + 1;
            
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
            rmaxC = max(rC(C > 0));
        end
        
        allPatterns{ri, ci} = C;
        allA{ri, ci} = A;
        allB{ri, ci} = B;
        fprintf('  Completed at t=%d\n', t);
    end
end

% === CREATE PUBLICATION FIGURE ===
% Layout parameters (normalized)
marginL = 0.08;
marginR = 0.01;
marginT = 0.08;
marginB = 0.01;
gapH = 0.01;
gapV = 0.015;

% Figure width in inches
figW = 8;

% Compute panel width (normalized)
panelW = (1 - marginL - marginR - (nCols-1)*gapH) / nCols;

% Compute available vertical space (normalized)
availableV = 1 - marginT - marginB - (nRows-1)*gapV;
panelH = availableV / nRows;

% Compute figure height so panels are square in inches
% panelW * figW = panelH * figH  =>  figH = panelW * figW / panelH
figH = panelW * figW / panelH;

figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 figW figH]);

for ri = 1:nRows
    for ci = 1:nCols
        C = allPatterns{ri, ci};
        A = allA{ri, ci};
        B = allB{ri, ci};
        
        left = marginL + (ci-1)*(panelW + gapH);
        bottom = 1 - marginT - ri*panelH - (ri-1)*gapV;
        
        ax = axes('Position', [left, bottom, panelW, panelH]);
        
        minC = min(C(C>0)); 
        maxC = max(C(C>0));
        dum1 = zeros(size(C));
        dum1(C>0) = 1 - (C(C>0)-minC) / (maxC-minC);
        img1 = cat(3, 0.8*A+dum1, 0.8*A+dum1*0.6+0.6*B, 0.8*A+B);
        
        imagesc(img1);
        axis equal off;

        % Column labels (top row only)
        if ri == 1
            title(sprintf('\\beta = %.2f', BETA_vals(ci)), ...
                'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial');
        end
        
        % Row labels (left column only)
        if ci == 1
            text(-0.05, 0.5, sprintf('k\\beta = %.1f', KB_vals(ri)), ...
                'Units', 'normalized', ...
                'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
                'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        end
    end
end

% Export
print(gcf, 'Fig_CA_kbeta_sweep.png', '-dpng', '-r300');
print(gcf, 'Fig_CA_kbeta_sweep.pdf', '-dpdf', '-vector');
fprintf('Saved: Fig_CA_kbeta_sweep.png and .pdf\n');