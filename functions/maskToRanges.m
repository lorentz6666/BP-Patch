function ranges = maskToRanges(mask)
mask = logical(mask(:));
d = diff([0; mask; 0]);
st = find(d==1);
ed = find(d==-1)-1;
ranges = [st ed];
end

