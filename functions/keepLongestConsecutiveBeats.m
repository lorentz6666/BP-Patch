function keep_beats2 = keepLongestConsecutiveBeats(keep_beats, min_consecutive_beats)
keep_beats = logical(keep_beats(:));
keep_beats2 = false(size(keep_beats));

if ~any(keep_beats)
    return;
end

d = diff([false; keep_beats; false]);
st = find(d == 1);
ed = find(d == -1) - 1;

run_len = ed - st + 1;
[max_len, kmax] = max(run_len);

if isempty(max_len) || max_len < min_consecutive_beats
    return;
end

keep_beats2(st(kmax):ed(kmax)) = true;
end

