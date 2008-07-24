% period: "Detect f=1/3"
% Tests for the presence of a sinusoid with frequency 1/3
% Usage:  plot_period(signal,display)
% (display=1): plot periodogram if sinusoid is detected

function [pvalue] = plot_period(signal, display)
    config;
    if ischar(signal),
        signal = get_signal(signal);
    end

    % Compute periodogram
    Y = fft(signal);
    N = length(Y);
    power = abs(Y).^2/N;

    % What's so special about the first value anyway?
    first = power(1);
    power = power(2:end);

    % We're testing to see if the points at 1/3 Hz and 2/3 Hz are
    % different from the points everywhere else. Indices are not in Hz
    % but in Hz * N, so multiply by N too.
    est = mean(power([N/3, 2*N/3]));
    F = (N-3)*est/(dot(signal,signal) - first - (2*est));
    pvalue = 1 - fcdf(F, 3-1, N-3);
    fprintf('p-value: %g\n\n', pvalue);

    if display,
        % 1+floor(N/2) because the remaining values are redundant. 2 instead
        % of 1 because a huge spike occurs at power(1) for no reason.
        % What's weird: we use power(1) in the above calculation.
        ind = 1:floor(N/2);
        plot(ind/N, power(ind));
        grid on;
        xlabel('Cycles per base')
        title('Power spectrum(FFT(free energy signal))')
    end
end
