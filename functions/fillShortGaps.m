function mask2 = fillShortGaps(mask, maxGap)
mask = logical(mask(:));
mask2 = mask;
d = diff([0; mask; 0]);
st0 = find(d==-1);
ed0 = find(d==1)-1;
gap_len = ed0 - st0 + 1;
for k = 1:numel(gap_len)
    if gap_len(k) <= maxGap
        left_on  = (st0(k) > 1) && mask(max(st0(k)-1,1));
        right_on = (ed0(k) < numel(mask)) && mask(min(ed0(k)+1,numel(mask)));
        if left_on && right_on
            mask2(st0(k):ed0(k)) = true;
        end
    end
end
end

