function [keep_seg, keep_beats] = simplePulseQC(amp_vals, jump_thresh_ratio, rise_ratio_thresh, min_beats_per_seg, peak_vals, trough_vals)
amp_vals = amp_vals(:);
peak_vals = peak_vals(:);
trough_vals = trough_vals(:);

n = numel(amp_vals);
keep_beats = true(n,1);

if n < min_beats_per_seg
    keep_seg = false;
    keep_beats = false(n,1);
    return;
end

med_amp = median(abs(amp_vals));
if med_amp <= 0 || ~isfinite(med_amp)
    keep_seg = false;
    keep_beats = false(n,1);
    return;
end

jump_ratio = [1; abs(diff(amp_vals)) / max(med_amp, eps)];
keep_beats(jump_ratio > jump_thresh_ratio) = false;

rise_ratio = peak_vals ./ max(peak_vals - trough_vals, eps);
keep_beats(rise_ratio < rise_ratio_thresh) = false;

if nnz(keep_beats) < min_beats_per_seg
    keep_seg = false;
    keep_beats(:) = false;
    return;
end

keep_seg = true;
end

