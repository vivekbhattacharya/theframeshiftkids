% I, given a file, call Kidnap for you.
%
% USAGE:
%     [signal, sequence] = get_signal('gene.txt')
function [signal, s] = get_signal(f)
global Config;

file = which(f);
if isempty(file)
    if exist(f) == 2, file = f;
    elseif strfind(f, 'http://'), file = f;
    else
        error(sprintf('get_signal cannot find "%s" in the path. Type which(''%s'') to confirm.', f, f));
    end
end

s = pearl('getseq.pl', ['"' file '"']);

signal_arg = '';
if Config.should_cache
    signal_arg = sprintf('auuccuccacuag "%s"', file);
else
    signal_arg = sprintf('-n auuccuccacuag "%s"', file);
end

% Now to tag on the temperature and choice of values
kelvin_temp = Config.temp + 273.15 % To account for Hao's idiotic idea to ask for the temperature in Kelvin instead of Celsius
signal_arg = ['-t ' num2str(kelvin_temp) ' -p Kidnap::' Config.values ' ' signal_arg];

output = pearl('scan_brightly.pl', signal_arg);

% Simulate load() on a string instead of a file.
signal = sscanf(output, '%f');
if isempty(signal)
    ensure = sprintf('I cannot pull signals. Ensure `perl -v` shows Perl is of version 5.6 or higher.\n');
    ensure = [ensure sprintf('Remember to run cleanpath(''path/to/BWFtools''); savepath.\n')];
    ensure = [ensure sprintf('Perl said for sequence: %s\n', s)];
    ensure = [ensure sprintf('Perl said for signal: %s\n', output)];
    error(ensure);
end

signal = [Config.signal_shift signal'];
