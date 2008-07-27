% Plot power spectrum for fft(free energy signal). Pass a non-zero
% value to display if you want a graph and p-value displayed.
%
% Usage:
% periodogram(signal, display)
%
% Returns the p-value from a most likely incorrectly calculated
% F-statistic. Use moving_periodogram to look at the p-values for a
% moving window of free energy signal values.
function [pvalue, snr_db] = periodogram(signal, display)
    config;
    if ischar(signal)
        signal = get_signal(signal);
    end
    
    % To make the length a multiple of 3
    % signal = signal(1:end-rem(length(signal),3));

    % Compute periodogram
    Y = fft(signal);
    N = length(Y);
    power = abs(Y).^2/N;

    % What's so special about the first value anyway? It represents the
    % power of the signal not having a frequency, which is
    % extraordinary sometimes.
    first = power(1);
    power = power(2:end);

    % We're testing to see if the points at 1/3 Hz and 2/3 Hz are
    % different from the points everywhere else. Indices are not in Hz
    % but in Hz * N, so multiply by N too.
    est = mean(power([N/3, 2*N/3]));
    F = (N-3)*est/(dot(signal,signal) - first - (2*est));
    pvalue = 1 - fcdf(F, 3-1, N-3);

    % Taken from est_par.
    x = ones(N, 1);
    y = (0:N-1)';
    regressor = [x sin(2*pi*y/3) cos(2*pi*y/3)];
    [b, b_int, r] = regress(signal', regressor);
    amp = sqrt((b(2)^2) + (b(3)^2));
    snr_db = 10*log10((amp^2/2) / var(r));

    if display
        fprintf('p-value: %g\n', pvalue);
        fprintf('SNR (dB): %g\n', snr_db);
        fprintf('\n');
    end

    if display
        ind = 1:floor(N/2);
        plot(ind/N, power(ind));
        grid on;
        xlabel('Cycles per base')
        title('Power spectrum(FFT(free energy signal))')
    end
end
