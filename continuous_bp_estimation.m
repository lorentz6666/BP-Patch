
clc;
clear;
close all;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'functions'));

path_data = "subject_pressure_waveform.csv";
gold_path = "reference_blood_pressure.csv";

gold = readtable(gold_path);
trim_n = 1;
gold(1:trim_n, :) = [];
gold.time_min = gold.time_min - gold.time_min(1);
cutoff_min = 100;
cutoff_left = 0;

SBPref = 113;
DBPref = 36;
alpha_map = 1/3;

data_airbag = readmatrix(path_data);

crop_start = 1;
data_airbag = data_airbag(crop_start:end,:);

input_column = 1;
if input_column < 1 || input_column > size(data_airbag, 2)
    error('input_column must be between 1 and %d.', size(data_airbag, 2));
end
CDC = data_airbag(:, input_column);
CDC = CDC(:);

Fs = 40;
N  = numel(CDC);

band = [20 60];
smooth_band_sec = 0.6;
smooth_env_sec  = 1.0;
dP_thr = -0.05;
minLen_sec = 5.0;
fillGap_sec = 0.25;
pad_sec = 0;

minPeakProm_all = 0.20;
minPeakDist_all = 15;

enable_simple_qc      = true;
jump_thresh_ratio     = 1.8;
rise_ratio_thresh     = 0.3;
min_beats_per_seg     = 4;
min_consecutive_beats = 8;


c_mode = getenv('BP_C_MODE');
if isempty(c_mode)
    c_mode = 'best_c';
end
use_best_c = strcmpi(c_mode, 'per_segment');
c_jump_ratio_thr = 1.5;
R2_try_value = 0.94;
wB = max(3, round(smooth_band_sec*Fs));
Pband = movmedian(CDC, wB);
Pband = movmean(Pband, wB);

wE = max(3, round(smooth_env_sec*Fs));
Penv = movmedian(CDC, wE);
Penv = movmean(Penv, wE);

dPenv = [0; diff(Penv)] * Fs;

mask_band    = (Pband >= band(1)) & (Pband <= band(2));
mask_deflate = (dPenv < dP_thr);
mask = mask_band & mask_deflate;

mask = fillShortGaps(mask, round(fillGap_sec*Fs));
ranges = maskToRanges(mask);

minLen = round(minLen_sec*Fs);
if ~isempty(ranges)
    ranges = ranges((ranges(:,2)-ranges(:,1)+1) >= minLen, :);
end

if isempty(ranges)
    error('自动选段 ranges 为空：请检查 band / dP_thr / 平滑窗口。');
end

ranges_idx = ranges;
pad = round(pad_sec * Fs);
ranges_idx(:,1) = max(1, ranges_idx(:,1) + pad);
ranges_idx(:,2) = min(N, ranges_idx(:,2) - pad);
valid_range = (ranges_idx(:,2) > ranges_idx(:,1));
ranges_idx = ranges_idx(valid_range,:);

if isempty(ranges_idx)
    error('加 pad 后所有选段都无效。');
end

numSeg = size(ranges_idx,1);

[~, pk_all] = findpeaks(CDC, 'MinPeakProminence', minPeakProm_all, 'MinPeakDistance', minPeakDist_all);

mirror = -CDC;
val_all = zeros(1, max(numel(pk_all)-1,0));
for k = 1:numel(pk_all)-1
    a = pk_all(k) + minPeakDist_all;
    b = pk_all(k+1);
    if a < b
        [~, mx] = max(mirror(a:b));
        val_all(k) = a + mx - 1;
    end
end
val_all = val_all(val_all > 0);

f_y              = cell(numSeg,1);
awave            = cell(numSeg,1);
temp_x           = cell(numSeg,1);
temp_mean        = cell(numSeg,1);
indexPulsePeak   = cell(numSeg,1);
indexPulseValley = cell(numSeg,1);
Peaks            = cell(numSeg,1);
Troughs          = cell(numSeg,1);

for i = 1:numSeg
    segStart = ranges_idx(i,1);
    segEnd   = ranges_idx(i,2);

    if segEnd - segStart < round(2*Fs)
        continue;
    end

    pk_idx = pk_all(pk_all>=segStart & pk_all<=segEnd);
    tr_idx = val_all(val_all>=segStart & val_all<=segEnd);

    if numel(tr_idx)<2 || numel(pk_idx)<2
        continue;
    end

    tx = floor((tr_idx(1:end-1) + tr_idx(2:end))/2);
    tm = arrayfun(@(j) mean(CDC(tr_idx(j):tr_idx(j+1))), 1:numel(tr_idx)-1).';

    if numel(tx) < 2
        continue;
    end

    tseg   = (tx(1):tx(end)).';
    fy_seg = interp1(tx, tm, tseg, 'spline', 'extrap');
    f_y{i} = [zeros(tx(1)-1,1); fy_seg; zeros(N - tx(end),1)];

    awave{i} = CDC - f_y{i};

    pk_use = pk_idx(pk_idx>=tseg(1) & pk_idx<=tseg(end));
    tr_use = tr_idx(tr_idx>=tseg(1) & tr_idx<=tseg(end));

    nPair = min(numel(pk_use), numel(tr_use));
    if nPair == 0
        continue;
    end

    pk_use = pk_use(1:nPair);
    tr_use = tr_use(1:nPair);

    peak_vals   = awave{i}(pk_use);
    trough_vals = awave{i}(tr_use);
    amp_vals    = peak_vals - trough_vals;

    if enable_simple_qc
        [keep_seg, keep_beats] = simplePulseQC(amp_vals, jump_thresh_ratio, rise_ratio_thresh, min_beats_per_seg, peak_vals, trough_vals);

        if ~keep_seg
            f_y{i}              = [];
            awave{i}            = [];
            temp_x{i}           = [];
            temp_mean{i}        = [];
            indexPulsePeak{i}   = [];
            indexPulseValley{i} = [];
            Peaks{i}            = [];
            Troughs{i}          = [];
            continue;
        end

        keep_beats = keepLongestConsecutiveBeats(keep_beats, min_consecutive_beats);

        if nnz(keep_beats) < min_consecutive_beats
            f_y{i}              = [];
            awave{i}            = [];
            temp_x{i}           = [];
            temp_mean{i}        = [];
            indexPulsePeak{i}   = [];
            indexPulseValley{i} = [];
            Peaks{i}            = [];
            Troughs{i}          = [];
            continue;
        end

        pk_use      = pk_use(keep_beats);
        tr_use      = tr_use(keep_beats);
        peak_vals   = peak_vals(keep_beats);
        trough_vals = trough_vals(keep_beats);
    end

    if numel(pk_use) < min_consecutive_beats || numel(tr_use) < min_consecutive_beats
        f_y{i}              = [];
        awave{i}            = [];
        temp_x{i}           = [];
        temp_mean{i}        = [];
        indexPulsePeak{i}   = [];
        indexPulseValley{i} = [];
        Peaks{i}            = [];
        Troughs{i}          = [];
        continue;
    end

    indexPulsePeak{i}   = pk_use(:).';
    indexPulseValley{i} = tr_use(:).';
    Peaks{i}            = peak_vals(:).';
    Troughs{i}          = trough_vals(:).';

    temp_x{i}    = tx(:);
    temp_mean{i} = tm(:);
end

Pc_vec = [];
Peak_vec = [];
Trough_vec = [];
PP_vec = [];
SegID_vec = [];
is_calib = false(0,1);

for i = 1:numSeg
    if isempty(f_y{i}) || isempty(awave{i}) || isempty(indexPulsePeak{i}) || isempty(indexPulseValley{i})
        continue;
    end

    pk_idx = indexPulsePeak{i};
    tr_idx = indexPulseValley{i};

    n = min([numel(pk_idx), numel(tr_idx)-1, numel(temp_mean{i})]);
    if n <= 0
        continue;
    end

    pk_idx = pk_idx(1:n);
    tr_idx = tr_idx(2:n+1);
    Pc_this = temp_mean{i}(1:n);

    Peak_this = awave{i}(pk_idx);
    Trough_this = abs(awave{i}(tr_idx));
    PP_this = Peak_this + Trough_this;
    Beat_this = pk_idx(:);

    valid_this = isfinite(Pc_this) & isfinite(Peak_this) & isfinite(Trough_this) & isfinite(PP_this) & ...
        (Peak_this > 0) & (Trough_this > 0) & (PP_this > 0);

    Pc_this     = Pc_this(valid_this);
    Peak_this   = Peak_this(valid_this);
    Trough_this = Trough_this(valid_this);
    PP_this     = PP_this(valid_this);
    Beat_this   = Beat_this(valid_this);

    if isempty(Pc_this)
        continue;
    end

    Pc_vec        = [Pc_vec; Pc_this(:)];
    Peak_vec      = [Peak_vec; Peak_this(:)];
    Trough_vec    = [Trough_vec; Trough_this(:)];
    PP_vec        = [PP_vec; PP_this(:)];
    SegID_vec     = [SegID_vec; repmat(i, numel(Pc_this), 1)];
    is_calib      = [is_calib; repmat(i==1, numel(Pc_this), 1)];
end

if isempty(Pc_vec)
    error('没有拼接出任何有效拍。请检查 ranges / 峰谷阈值 / 连续波数要求。');
end

best_seg_cal = NaN;
best_c = NaN;
best_phi = NaN;
best_R2 = -Inf;
best_Pc_cal = [];
best_PP_cal = [];
best_x_fit = [];
best_y_fit = [];
best_coeff = [];
c_try_all = NaN(numSeg,1);
R2_try_all = NaN(numSeg,1);

for s = 1:numSeg
    idxs = (SegID_vec == s);
    Pc_try = Pc_vec(idxs);
    PP_try = PP_vec(idxs);

    valid_try = isfinite(Pc_try) & isfinite(PP_try) & (PP_try > 0);
    Pc_try = Pc_try(valid_try);
    PP_try = PP_try(valid_try);

    if numel(PP_try) < 6
        continue;
    end

    x_try = Pc_try(3:end);
    y_try = log(PP_try(3:end));

    if numel(x_try) < 4
        continue;
    end

    coeff_try = polyfit(x_try, y_try, 1);
    c_try = coeff_try(1);
    intercept_try = coeff_try(2);

    if ~isfinite(c_try) || c_try <= 0
        continue;
    end

    term_inside = exp(-c_try * DBPref) - exp(-c_try * SBPref);
    if term_inside <= 0 || ~isfinite(term_inside)
        continue;
    end

    phi_try = intercept_try - log(term_inside);

    yhat_try = polyval(coeff_try, x_try);
    SS_res = sum((y_try - yhat_try).^2);
    SS_tot = sum((y_try - mean(y_try)).^2);
    R2_try = 1 - SS_res / max(SS_tot, eps);
    R2_try_all(s) = R2_try;
    if R2_try > R2_try_value
        c_try_all(s) = c_try;
    end
    if isfinite(R2_try) && (R2_try > best_R2)
        best_seg_cal = s;
        best_c = c_try;
        best_phi = phi_try;
        best_R2 = R2_try;
        best_Pc_cal = Pc_try;
        best_PP_cal = PP_try;
        best_x_fit = x_try;
        best_y_fit = y_try;
        best_coeff = coeff_try;
    end
end

if isnan(best_seg_cal)
    error('所有分段都无法得到正的标定斜率 c。请检查 band 或信号方向。');
end

c = best_c;
phi = best_phi;
R2_ref = best_R2;
Pc_cal = best_Pc_cal;
PP_cal = best_PP_cal;
x_fit = best_x_fit;
y_fit = best_y_fit;
coeff = best_coeff;

if use_best_c
    c_use_seg = repmat(best_c, numSeg, 1);
else
    c_use_seg = NaN(numSeg, 1);
    last_valid_c = best_c;

    for s = 1:numSeg
        c_candidate = c_try_all(s);
        is_bad_c = ~isfinite(c_candidate) || c_candidate <= 0;

        if ~is_bad_c && isfinite(last_valid_c) && last_valid_c > 0
            if (c_candidate > last_valid_c * c_jump_ratio_thr) || (c_candidate < last_valid_c / c_jump_ratio_thr)
                is_bad_c = true;
            end
        end

        if is_bad_c
            c_use_seg(s) = last_valid_c;
        else
            c_use_seg(s) = c_candidate;
            last_valid_c = c_candidate;
        end
    end
end

sbp_est = nan(numSeg,1);
dbp_est = nan(numSeg,1);
map_est = nan(numSeg,1);
mean_pp = nan(numSeg,1);
objval  = nan(numSeg,1);
nbeats  = nan(numSeg,1);
pp_cv   = nan(numSeg,1);
pc_std  = nan(numSeg,1);
bad_flag = false(numSeg,1);

for i = 1:numSeg
    idx = (SegID_vec == i);

    if nnz(idx) < min_consecutive_beats
        bad_flag(i) = true;
        continue;
    end

    Pc_i = Pc_vec(idx);
    Peak_i = Peak_vec(idx);
    Trough_i = Trough_vec(idx);
    PP_i = PP_vec(idx);
    c_i = c_use_seg(i);

    if ~isfinite(c_i) || c_i <= 0
        bad_flag(i) = true;
        continue;
    end

    ae_list = log(Peak_i)   - phi - c_i * Pc_i;
    be_list = log(PP_i)     - phi - c_i * Pc_i;
    ce_list = log(Trough_i) - phi - c_i * Pc_i;

    valid = isfinite(ae_list) & isfinite(be_list) & isfinite(ce_list);
    ae_list = ae_list(valid);
    be_list = be_list(valid);
    ce_list = ce_list(valid);
    PP_i = PP_i(valid);
    Pc_i = Pc_i(valid);

    if numel(ae_list) < min_consecutive_beats
        bad_flag(i) = true;
        continue;
    end

    x0 = [SBPref; DBPref];
    lb = [0; 0];
    ub = [220; 130];

    objfun = @(x) total_residual_local(x, ae_list, be_list, ce_list, c_i, alpha_map, SBPref, DBPref);

    if license('test','optimization_toolbox')
        options = optimoptions('fmincon', 'Display','off', ...
            'TolFun',1e-10, 'MaxFunctionEvaluations',5000, 'Algorithm','interior-point');
        [x, fval] = fmincon(objfun, x0, [], [], [], [], lb, ub, [], options);
    else
        penalty_fun = @(x) objfun(project_box_local(x,lb,ub));
        x = fminsearch(penalty_fun, x0, optimset('Display','off','MaxFunEvals',5000,'MaxIter',3000));
        x = project_box_local(x,lb,ub);
        fval = penalty_fun(x);
    end

    SBP = x(1);
    DBP = x(2);
    MAP = DBP + alpha_map*(SBP-DBP);

    sbp_est(i) = SBP;
    dbp_est(i) = DBP;
    map_est(i) = MAP;
    mean_pp(i) = mean(PP_i);
    objval(i)  = fval;
    nbeats(i)  = numel(ae_list);
    pp_cv(i)   = std(PP_i) / max(mean(PP_i), eps);
    pc_std(i)  = std(Pc_i);


end

first_sbp_idx = find(isfinite(sbp_est), 1, 'first');
first_dbp_idx = find(isfinite(dbp_est), 1, 'first');

if isempty(first_sbp_idx) || isempty(first_dbp_idx)
    error('sbp_est 或 dbp_est 全部为 NaN，无法做参考值平移。');
end

sbp_est = sbp_est + SBPref - sbp_est(first_sbp_idx);
dbp_est = dbp_est + DBPref - dbp_est(first_dbp_idx);

map_est = dbp_est + alpha_map*(sbp_est-dbp_est);

finalBoxW_cm_target = 7.48;
finalBoxH_cm_target = 2.75;
finalFsLabel_target = 10;
finalFsTick_target  = 6;
finalAxLW_target    = 0.55;
finalTickLen_cm     = 0.05;

col_sbp = [241 64 64] / 255;
col_map = [167 206 115] / 255;
col_dbp = [26 111 223] / 255;

fsLabel = finalFsLabel_target;
fsTick  = finalFsTick_target;
axLW    = finalAxLW_target;
lineLW  = 1.00;
goldLW  = 0.75;

bpBoxW_cm = finalBoxW_cm_target;
bpBoxH_cm = finalBoxH_cm_target;
tickLenNorm = finalTickLen_cm / max(bpBoxW_cm, bpBoxH_cm);

leftMargin_cm   = 1.05;
bottomMargin_cm = 0.75;
rightMargin_cm  = 0.18;
topMargin_cm    = 0.12;

figW_cm = leftMargin_cm + bpBoxW_cm + rightMargin_cm;
figH_cm = bottomMargin_cm + bpBoxH_cm + topMargin_cm;

fig = figure('Name','Estimated BP summary0','Color','w', ...
    'Units','centimeters', ...
    'Position',[5 5 figW_cm figH_cm], ...
    'PaperUnits','centimeters', ...
    'PaperPosition',[0 0 figW_cm figH_cm], ...
    'PaperSize',[figW_cm figH_cm], ...
    'Renderer','painters');

gold0 = readtable(gold_path);
gold0(1:trim_n, :) = [];

ax1 = axes(fig, ...
    'Units','centimeters', ...
    'Position',[leftMargin_cm bottomMargin_cm bpBoxW_cm bpBoxH_cm]);
hold(ax1, 'on');

seg_valid = ~isnan(sbp_est) & ~isnan(dbp_est) & ~isnan(map_est);
seg_time_min = ((ranges_idx(:,1) + ranges_idx(:,2)) / 2) / Fs / 60;
time_keep = (seg_time_min >= cutoff_left) & (seg_time_min <= cutoff_min);

x = seg_time_min(seg_valid & time_keep);
if numel(x) < 2
    error('Valid BP points are fewer than 2; pchip interpolation cannot be performed.');
end

x = x - cutoff_left;
y_sbp = sbp_est(seg_valid & time_keep);
y_dbp = dbp_est(seg_valid & time_keep);
y_map = map_est(seg_valid & time_keep);

[x, sort_idx] = sort(x(:));
y_sbp = y_sbp(sort_idx);
y_dbp = y_dbp(sort_idx);
y_map = y_map(sort_idx);

[x, unique_idx] = unique(x, 'stable');
y_sbp = y_sbp(unique_idx);
y_dbp = y_dbp(unique_idx);
y_map = y_map(unique_idx);

if numel(x) < 2
    error('Valid BP points are fewer than 2 after de-duplication; pchip interpolation cannot be performed.');
end

x_smooth = linspace(min(x), max(x), 500);
sbp_smooth = interp1(x, y_sbp, x_smooth, 'pchip');
dbp_smooth = interp1(x, y_dbp, x_smooth, 'pchip');
map_smooth = interp1(x, y_map, x_smooth, 'pchip');

gold0.time_min = gold0.time_min - gold0.time_min(1);
gold_keep = (gold0.time_min >= cutoff_left) & (gold0.time_min <= cutoff_min);
gold0 = gold0(gold_keep, :);
gold0.time_min = gold0.time_min - cutoff_left;

if isempty(gold0)
    common_end_min = max(x_smooth);
    x_plot_keep = x_smooth <= common_end_min;
    gold_plot = gold0;
else
    common_end_min = min(max(x_smooth), max(gold0.time_min));
    if ~isfinite(common_end_min) || common_end_min <= 0
        common_end_min = max(x_smooth);
    end
    x_plot_keep = x_smooth <= common_end_min;
    gold_plot = gold0(gold0.time_min <= common_end_min, :);
end

plot(ax1, x_smooth(x_plot_keep), sbp_smooth(x_plot_keep), '-', ...
    'Color', col_sbp, 'LineWidth', lineLW, 'HandleVisibility', 'off');
plot(ax1, x_smooth(x_plot_keep), map_smooth(x_plot_keep), '-', ...
    'Color', col_map, 'LineWidth', lineLW, 'HandleVisibility', 'off');
plot(ax1, x_smooth(x_plot_keep), dbp_smooth(x_plot_keep), '-', ...
    'Color', col_dbp, 'LineWidth', lineLW, 'HandleVisibility', 'off');

if ~isempty(gold_plot)
    plot(ax1, gold_plot.time_min, gold_plot.SBP, '-', ...
        'Color', [0.70 0.70 0.70], 'LineWidth', goldLW, 'HandleVisibility', 'off');
    p2 = plot(ax1, gold_plot.time_min, gold_plot.DBP, '-', ...
        'Color', [0.70 0.70 0.70], 'LineWidth', goldLW);
    p2.HandleVisibility = 'off';
    p3 = plot(ax1, gold_plot.time_min, gold_plot.MAP, '-', ...
        'Color', [0.70 0.70 0.70], 'LineWidth', goldLW);
    p3.HandleVisibility = 'off';
end

xlabel(ax1, 'Time (min)', ...
    'FontName', 'Arial', ...
    'FontSize', fsLabel, ...
    'FontWeight', 'normal');
ylabel(ax1, 'BP (mmHg)', ...
    'FontName', 'Arial', ...
    'FontSize', fsLabel, ...
    'FontWeight', 'normal');

yDataForLim = [sbp_smooth(x_plot_keep), map_smooth(x_plot_keep), dbp_smooth(x_plot_keep)];
if ~isempty(gold_plot)
    yDataForLim = [yDataForLim, gold_plot.SBP(:).', gold_plot.MAP(:).', gold_plot.DBP(:).'];
end
yDataForLim = yDataForLim(isfinite(yDataForLim));

if isempty(yDataForLim)
    yLim = [50 170];
else
    yMin = min(yDataForLim);
    yMax = max(yDataForLim);
    if yMin == yMax
        yMin = yMin - 10;
        yMax = yMax + 10;
    end
    yPad = max(5, 0.08 * (yMax - yMin));
    yLim = [floor((yMin - yPad) / 10) * 10, ceil((yMax + yPad) / 10) * 10];
    if diff(yLim) < 20
        yMid = mean(yLim);
        yLim = [floor((yMid - 10) / 10) * 10, ceil((yMid + 10) / 10) * 10];
    end
end

ylim(ax1, yLim);
xlim(ax1, [0 common_end_min]);

if common_end_min <= 10
    xTicks = linspace(0, common_end_min, 3);
else
    xTickStep = 5;
    xTicks = 0:xTickStep:(floor(common_end_min / xTickStep) * xTickStep);
    if numel(xTicks) < 3
        xTicks = linspace(0, common_end_min, 3);
    end
end
xTicks = unique(round(xTicks));
xTicks = xTicks(xTicks >= 0 & xTicks <= common_end_min);
xticks(ax1, xTicks);
xtickformat(ax1, '%.0f');

yTickStep = max(10, ceil(diff(yLim) / 3 / 10) * 10);
yTicks = ceil(yLim(1) / yTickStep) * yTickStep:yTickStep:floor(yLim(2) / yTickStep) * yTickStep;
if numel(yTicks) < 3
    yTicks = linspace(yLim(1), yLim(2), 3);
end
yTicks = unique(round(yTicks));
yticks(ax1, yTicks);
ytickformat(ax1, '%.0f');

set(ax1, ...
    'FontName', 'Arial', ...
    'FontSize', fsTick, ...
    'Box', 'off', ...
    'LineWidth', axLW, ...
    'TickDir', 'out', ...
    'TickLength', [tickLenNorm tickLenNorm], ...
    'Layer', 'top', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'Color', 'none');
ax1.XAxis.Exponent = 0;

pos = get(ax1, 'Position');
ax2 = axes(fig, ...
    'Units','centimeters', ...
    'Position', pos, ...
    'Color', 'none', ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'XLim', ax1.XLim, ...
    'YLim', ax1.YLim, ...
    'XTick', [], ...
    'YTick', [], ...
    'Box', 'on', ...
    'LineWidth', axLW, ...
    'TickLength', [tickLenNorm tickLenNorm], ...
    'Layer', 'top', ...
    'HitTest', 'off');
ax2.XColor = ax1.XColor;
ax2.YColor = ax1.YColor;

linkaxes([ax1 ax2], 'xy');
axes(ax1);

set(fig, ...
    'Units', 'centimeters', ...
    'Position', [5 5 figW_cm figH_cm], ...
    'PaperUnits', 'centimeters', ...
    'PaperPosition', [0 0 figW_cm figH_cm], ...
    'PaperSize', [figW_cm figH_cm]);

outFolder = 'Figures_PNG';
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

thisMFile = mfilename;
if isempty(thisMFile)
    thisMFile = 'BP_summary';
end

outName = fullfile(outFolder, [thisMFile '.png']);
exportgraphics(fig, outName, ...
    'Resolution', 600, ...
    'BackgroundColor', 'white');



seg_time_min_shift = seg_time_min - cutoff_left;
gold_sbp = interp1(gold.time_min, gold.SBP, seg_time_min_shift, 'linear', 'extrap');
gold_dbp = interp1(gold.time_min, gold.DBP, seg_time_min_shift, 'linear', 'extrap');
gold_map = interp1(gold.time_min, gold.MAP, seg_time_min_shift, 'linear', 'extrap');

idx_eval = seg_valid & time_keep & isfinite(gold_sbp) & isfinite(gold_dbp) & isfinite(gold_map);

err_sbp = sbp_est(idx_eval) - gold_sbp(idx_eval);
err_dbp = dbp_est(idx_eval) - gold_dbp(idx_eval);
err_map = map_est(idx_eval) - gold_map(idx_eval);

fprintf('\n==== 与金标准对比 (按段中位时间) ====\n');
fprintf('SBP: mean=%.3f, var=%.3f, MAE=%.3f\n', mean(err_sbp), var(err_sbp), mean(abs(err_sbp)));
fprintf('DBP: mean=%.3f, var=%.3f, MAE=%.3f\n', mean(err_dbp), var(err_dbp), mean(abs(err_dbp)));
fprintf('MAP: mean=%.3f, var=%.3f, MAE=%.3f\n', mean(err_map), var(err_map), mean(abs(err_map)));
