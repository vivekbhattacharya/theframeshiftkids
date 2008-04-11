% I, given a file, call Kidnap for you.
%
% USAGE:
%     [signal, sequence] = get_signal('gene.txt')
function [signal, s] = get_signal(f)

file = which(f);
if isempty(file)
    if exist(f) == 2, file = f;
    elseif strfind(f, 'http://'), file = f;
    else
        error(sprintf('get_signal cannot find "%s" in the path. Type which(''%s'') to confirm.', f, f));
    end
end

s = pearl('getseq.pl', ['"' file '"']);
signal = pearl('scan_brightly.pl', sprintf('auuccuccacuag "%s"', file));

% Simulate load() on a string instead of a file.
signal = sscanf(signal, '%f');
if isempty(signal)
    ensure = sprintf('I cannot pull signals. Ensure `perl -v` shows Perl is of version 5.6 or higher.\n');
    ensure = [ensure sprintf('Remember to run cleanpath(''path/to/BWFtools''); savepath.\n')];
    ensure = [ensure sprintf('Perl said for sequence: %s\n', s)];
    ensure = [ensure sprintf('Perl said for signal: %s\n', signal)];
    error(ensure);
end

global Config;
signal = [Config.signal_shift signal'];
