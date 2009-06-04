% I, given a file, call Kidnap for you. Returns a vector of free
% energies (dG) at each nucleotide during mRNA translation.
%
% USAGE:
%     [signal, sequence] = get_signal('gene.txt')
function [signal, s] = get_signal(f)
    config;
    global Config;

    file = superwhich(f);
    s = pearl('getseq.pl', ['"' file '"']);

    cache_arg = '-n';
    if Config.should_cache
        cache_arg = '';
    end

    args = sprintf('-t %g -p "%s" %s %s "%s"', ...
                   Config.temp, Config.energies, ...
                   cache_arg, Config.tail, file);
    output = pearl(Config.scanner, args);

    % Simulate load() on a string instead of a file.
    signal = sscanf(output, '%f');
    if isempty(signal)
        ensure = sprintf('I cannot pull signals. Ensure `perl -v` shows Perl is of version 5.6 or higher.\n');
        ensure = [ensure sprintf('Remember to run cleanpath(''path/to/BWFtools''); savepath.\n')];
        ensure = [ensure sprintf('Perl said for sequence: %s\n', s)];
        ensure = [ensure sprintf('Perl said for signal: %s\n', output)];
        error(ensure);
    end

    signal = signal';
    shift = Config.signal_shift;

    if shift > 0
        signal = [zeros(1, shift) signal];
    elseif shift < 0
        signal = [signal((1 - shift):end) zeros(1, -shift)];
    end
end
