function inspect_one_segment(CDC, range_seg, seg_id, minPeakProm, minPeakDist)
pulseWave = CDC(range_seg(1):range_seg(2));

[~, indexPulsePeak] = findpeaks(pulseWave, ...
    'MinPeakProminence', minPeakProm, ...
    'MinPeakDistance', minPeakDist);

mirror = -pulseWave;
indexPulseValley = zeros(1, max(numel(indexPulsePeak)-1,0));
for k = 1:numel(indexPulsePeak)-1
    a = indexPulsePeak(k) + minPeakDist;
    b = indexPulsePeak(k+1);
    if a < b
        [~, mx] = max(mirror(a:b));
        indexPulseValley(k) = a + mx - 1;
    end
end
indexPulseValley = indexPulseValley(indexPulseValley>0);

if numel(indexPulsePeak) < 3 || numel(indexPulseValley) < 3
    figure('Color','w');
    plot(pulseWave, 'k');
    title(sprintf('Segment %d: not enough peaks/valleys', seg_id));
    grid on;
    return;
end

temp_mean = [];
temp_x = [];
for j = 1:length(indexPulseValley)-1
    temp_mean(j,1) = mean(pulseWave(indexPulseValley(j):indexPulseValley(j+1),1));
    temp_x(j,1) = floor((indexPulseValley(j)+indexPulseValley(j+1))/2);
end

t = (temp_x(1):temp_x(end))';
f_y = interp1(temp_x, temp_mean, t, 'spline');
zeros_front = zeros(temp_x(1)-1, 1);
zeros_back = zeros(length(pulseWave)-temp_x(end), 1);
f_y = [zeros_front; f_y; zeros_back];

if length(pulseWave) ~= length(f_y)
    n = min(length(pulseWave), length(f_y));
    pulseWave = pulseWave(1:n);
    f_y = f_y(1:n);
end

awave = pulseWave - f_y;

tempPeakindex = (indexPulsePeak>=temp_x(1) & indexPulsePeak<=temp_x(end));
tempTroughindex = (indexPulseValley>=temp_x(1) & indexPulseValley<=temp_x(end));

Peaks = awave(indexPulsePeak(tempPeakindex));
Troughs = awave(indexPulseValley(tempTroughindex));

n = min([numel(Peaks), numel(Troughs)-1, numel(temp_mean)]);
if n > 0
    peak_vals = Peaks(1:n);
    trough_vals = abs(Troughs(2:n+1));
    PP_vals = peak_vals + trough_vals;
    Pc_vals = temp_mean(1:n);
else
    peak_vals = [];
    trough_vals = [];
    PP_vals = [];
    Pc_vals = [];
end

figure('Name', sprintf('Inspect segment %d', seg_id), 'Color', 'w');

subplot(3,1,1); hold on; grid on;
plot(pulseWave, 'k', 'DisplayName', 'Raw waveform');
plot(f_y, 'r-', 'LineWidth',1.5, 'DisplayName', 'Static fit');
plot(indexPulsePeak(tempPeakindex), pulseWave(indexPulsePeak(tempPeakindex)), ...
    'mo', 'MarkerFaceColor','m', 'DisplayName','Peaks');
plot(indexPulseValley(tempTroughindex), pulseWave(indexPulseValley(tempTroughindex)), ...
    'cv', 'MarkerFaceColor','c', 'DisplayName','Valleys');
plot(temp_x, temp_mean, 'bo', 'DisplayName', 'Mean points');
title(sprintf('Segment %d: raw + static fit', seg_id));
legend('Location','best');

subplot(3,1,2); hold on; grid on;
plot(awave, 'b', 'DisplayName', 'AC component');
plot(indexPulsePeak(tempPeakindex), awave(indexPulsePeak(tempPeakindex)), ...
    'ro', 'MarkerFaceColor','r', 'DisplayName','Peak AC');
plot(indexPulseValley(tempTroughindex), awave(indexPulseValley(tempTroughindex)), ...
    'go', 'MarkerFaceColor','g', 'DisplayName','Trough AC');
yline(0,'--k');
title('AC component');
legend('Location','best');

subplot(3,1,3); hold on; grid on;
if ~isempty(PP_vals)
    plot(Pc_vals, peak_vals, 'ro-', 'DisplayName','Peak');
    plot(Pc_vals, trough_vals, 'go-', 'DisplayName','Trough');
    plot(Pc_vals, PP_vals, 'bo-', 'DisplayName','PP');
end
xlabel('Pc');
ylabel('Amplitude');
title('Beat features in this segment');
legend('Location','best');

fprintf('\n---- Segment %d inspection ----\n', seg_id);
if ~isempty(PP_vals)
    fprintf('n beats used   = %d\n', numel(PP_vals));
    fprintf('mean PP        = %.4f\n', mean(PP_vals));
    fprintf('std PP         = %.4f\n', std(PP_vals));
    fprintf('PP CV          = %.4f\n', std(PP_vals)/max(mean(PP_vals),eps));
    fprintf('mean Pc        = %.4f\n', mean(Pc_vals));
    fprintf('std Pc         = %.4f\n', std(Pc_vals));
end
end
