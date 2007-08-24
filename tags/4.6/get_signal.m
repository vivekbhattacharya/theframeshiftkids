function [Signal, S] = get_signal(f)

file = which(f);
if isempty(file)
    if exist(f) == 2, file = f;
    elseif strfind(f, 'http://'), file = f;
    else
        error(sprintf('get_signal cannot find ``%s" in your path. Type which(''%s'') to confirm.', f, f));
    end;
end;

if(isempty('perl')); error(['I cannot find Perl 5.6 or above.' ...
   'Please download it or add it to your PATH (see website).']); end;

% Load prfb, and generate a fasta file with column width of 60.
S = getseq(file);

% Run free2bind. Make sure to include its directory into -I.
Include = fileparts(which('Smooth.pm'));
   Template = 'perl -I"%s" "%s" auuccuccacuag "%s"';
   Command = sprintf(Template, Include, which('scan_brightly.pl'), file);
[status, raw_signal] = system(Command);

% Simulate load() on a string instead of a file.
Signal = str2num(raw_signal);
Signal = Signal';
if(isempty(Signal))
    ensure = 'I cannot pull signals. Ensure `perl` is outputting the rite stuff.';
    ensure = [ensure sprintf('\n') 'Also, ensure Perl is of version 5.6 or above.'];
    ensure = [ensure sprintf('\n') sprintf('\n') 'Perl says: ' raw_signal];
    ensure = [ensure sprintf('\n') 'Your sequence was: ' S];
    error(ensure);
end