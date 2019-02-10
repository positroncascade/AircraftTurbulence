function smooth_spectrum = smoothing_filter(spectrum)
    smooth_spectrum = 0.25 * spectrum(1:end-2) + 0.5 * spectrum(2:end-1) + 0.25 * spectrum(3:end);
    smooth_spectrum = [spectrum(1); smooth_spectrum; spectrum(end)];
end

