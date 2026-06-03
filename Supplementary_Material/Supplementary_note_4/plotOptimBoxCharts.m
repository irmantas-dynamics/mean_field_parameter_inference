function plotOptimBoxCharts(dirnames, varargin)
% PLOTOPTIMBOXCHARTS  Plot boxcharts from optimization results across multiple directories.
%
% USAGE:
%   plotOptimBoxCharts(dirnames)
%   plotOptimBoxCharts(dirnames, 'ParamName', ParamValue, ...)
%
% REQUIRED INPUT:
%   dirnames   - cell array of directory paths, e.g. {'dir1', '../dir2'}
%
% OPTIONAL NAME-VALUE PAIRS:
%   'prm'        - true parameter vector         (default: [1;20;3.25;1;10])
%   'lbl'        - x-axis tick labels (cell)     (default: {'200','500','1000','10000','100000','inf'})
%   'lbl2'       - LaTeX x-axis tick labels      (default: {'2\cdot10^2',...,'\infty'})
%   'pv'         - LaTeX y-axis labels per param (default: delta, J, etab, beta, taum)
%   'filePattern'- sprintf pattern with one %s   (default: '/QIF_AD_%s_optim_results.csv')
%   'bxWidth'    - box width                     (default: 0.2)
%   'deviation'  - x-offset between groups       (default: 0.15)
%   'colors'     - nDirs x 3 color matrix        (default: red→blue ramp)
%   'legendNames'- cell array of legend labels   (default: {'Dir1','Dir2',...})
%   'ax1'        - subplot Position matrix       (default: built-in 5-panel layout)
%
% EXAMPLE (replicates original 2-directory script):
%   plotOptimBoxCharts({'', '../optimization_results'}, ...
%       'legendNames', {'PSO','DE'})

% ── parse inputs ──────────────────────────────────────────────────────────────
p = inputParser;
nDirs = numel(dirnames);

addRequired(p,  'dirnames');
addParameter(p, 'prm',         [1;20;3.25;1;10]);
addParameter(p, 'lbl',         {'200','500','1000','10000','100000','inf'});
addParameter(p, 'lbl2',        {'2\cdot 10^2','5\cdot 10^2','10^3','10^4','10^5','\infty'});
addParameter(p, 'pv',          {'$\delta\Delta$','$\delta J$','$\delta\bar{\eta}$', ...
                                 '$\delta\tau_m$','$\delta\tau_a$'});
addParameter(p, 'filePattern', '/QIF_AD_%s_optim_results.csv');
addParameter(p, 'bxWidth',     0.2);
addParameter(p, 'deviation',   0.15);
addParameter(p, 'colors',      defaultColors(nDirs));
addParameter(p, 'legendNames', arrayfun(@(k) sprintf('Dir%d',k), 1:nDirs, 'UniformOutput', false));
addParameter(p, 'ax1',         defaultAxPos());

parse(p, dirnames, varargin{:});
opt = p.Results;

nFiles  = numel(opt.lbl);
nParams = numel(opt.pv);

% x-offsets: centre the group of nDirs boxes around each tick
offsets = linspace(-(nDirs-1)/2, (nDirs-1)/2, nDirs) * opt.deviation;

% storage for legend handles
paramHandles = cell(nDirs, nParams);   % paramHandles{d,i}
lossHandles  = cell(nDirs, 1);         % lossHandles{d}

% ── main loop over directories ────────────────────────────────────────────────
for idir = 1:nDirs
    col  = opt.colors(idir,:);
    offs = offsets(idir);

    mnm = ones(nParams,1) *  1e4;
    mxm = ones(nParams,1) * -1e4;

    for ifile = 1:nFiles
        if contains(dirnames{idir},'Invasive')
            fname = sprintf([dirnames{idir}, opt.filePattern(1:end-4),'_ext.csv'], opt.lbl{ifile});
        else
            fname = sprintf([dirnames{idir}, opt.filePattern], opt.lbl{ifile});
        end
        
        M     = readmatrix(fname);
        tmp   = ones(size(M,1), 1);
        xPos  = ifile + offs;

        % ── parameter figure ──────────────────────────────────────────────
        figure(1);
        for i1 = 1:nParams
            subplot(3, 2, i1);
            tmpY = (M(:,i1) - opt.prm(i1)) / opt.prm(i1);

            h = boxchart(tmp*xPos, tmpY);
            hold on;
            h.BoxFaceColor = col;
            h.BoxWidth     = opt.bxWidth;
            h.MarkerStyle  = '.';
            h.MarkerColor  = col;
            ylabel(opt.pv{i1}, 'Interpreter', 'latex');

            if ifile == 1          % store first-file handle for legend
                paramHandles{idir,i1} = h;
            end

            mnm(i1) = min(mnm(i1), min(tmpY));
            mxm(i1) = max(mxm(i1), max(tmpY));
        end

        % ── loss figure ───────────────────────────────────────────────────
        figure(2);
        hL = boxchart(tmp*xPos, M(:,end) + 1e-7);
        hold on;
        hL.BoxFaceColor = col;
        hL.BoxWidth     = opt.bxWidth;
        hL.MarkerStyle  = '.';
        hL.MarkerColor  = col;
        xlabel('$N$',    'Interpreter', 'latex');
        ylabel('Loss',   'Interpreter', 'latex');
        hold on;
        box on;

        if ifile == 1
            lossHandles{idir} = hL;
        end
    end % ifile
end % idir

% ── cosmetic: parameter figure ────────────────────────────────────────────────
figure(1);
for i1 = 1:nParams
    subplot(3, 2, i1);
    plot([0.5, nFiles+0.5], [0, 0], ':', 'Color', [0.4 0.4 0.4]);
    hold on;
    xlim([0.5, nFiles+0.5]);
    xticks(1:nFiles+1);

    if i1 < 3
        set(gca, 'XTickLabel', []);
    else
        set(gca, 'XTickLabel', opt.lbl2);
    end

    legend([paramHandles{:,i1}], opt.legendNames, 'Interpreter', 'latex');
    box on;
end

if size(opt.ax1,1) >= nParams
    for i1 = 1:nParams
        subplot(3, 2, i1);
        set(gca, 'Position', opt.ax1(i1,:));
    end
end

set(gcf, 'PaperSize', [17 9]);
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [100, 1000, 17, 9]);

% ── cosmetic: loss figure ─────────────────────────────────────────────────────
figure(2);
set(gca, 'XTickLabel', opt.lbl2);
set(gca, 'YScale', 'log');
xlim([0.5, nFiles+0.5]);
xticks(1:nFiles+1);
legend([lossHandles{:}], opt.legendNames, 'Interpreter', 'latex');

set(gcf, 'PaperSize', [10 6]);
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [100, 1000, 10, 6]);

end % function


% ── helpers ───────────────────────────────────────────────────────────────────
function C = defaultColors(n)
% Interpolate from red to blue for n directories.
    if n == 1
        C = [0.9 0.25 0.25];
        return;
    end
    t = linspace(0, 1, n)';
    C = [0.9*(1-t) + 0.0*t, ...   % R
         0.25*(1-t) + 0.15*t, ...  % G
         0.25*(1-t) + 0.95*t];     % B
end

function ax1 = defaultAxPos()
    ax1 = [0.07, 0.72, 0.40, 0.26;
           0.56, 0.72, 0.40, 0.26;
           0.07, 0.44, 0.40, 0.26;
           0.56, 0.44, 0.40, 0.26;
           0.27, 0.10, 0.40, 0.26];
end
