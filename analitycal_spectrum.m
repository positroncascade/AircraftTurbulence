function spectrum = analitycal_spectrum(system, w, input_num)
    total_response = bode(system, w);
    turbulence_response = total_response(:, input_num, :);
    spectrum = turbulence_response .* turbulence_response;
end