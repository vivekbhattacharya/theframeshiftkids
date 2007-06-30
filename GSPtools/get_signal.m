function [Signal, S] = get_signal(f)

file = which(f);
if isempty(file)
    if exist(f) == 2, file = f;
    elseif strfind(f, 'http://'), file = f;
    else
        error(['Cannot find ``' file '" in your path']);
    end;
end;

if(isempty('perl.exe')); error(['I cannot find Perl 5.6 or above.' ...
   'Please download it to the VCL C: drive.']); end;

% Load prfb, and generate a fasta file with column width of 60.
S = getseq(file); Fasta = tempname;
write2fasta(Fasta, S, 'prfb', 60);

% Run free2bind. Make sure to include its directory into -I.
Include = fileparts(which('FreeAlign.pm'));
   Template = 'perl.exe -I"%s" "%s" -e -q -p FREIER auuccuccacuag "%s"';
   Command = sprintf(Template, Include, which('free_scan.pl'), Fasta);
[status, raw_signal] = system(Command);

% Simulate load() on a string instead of a file.
Signal = str2num(raw_signal);
Signal = Signal';
if(isempty(Signal))
    ensure = 'I cannot pull signals. Ensure `perl.exe` is outputting the rite stuff.';
    ensure = [ensure sprintf('\n') 'Also, ensure Perl is of version 5.6 or above.'];
    ensure = [ensure sprintf('\n') 'Also, ensure free2bind is in the path.'];
    ensure = [ensure sprintf('\n') sprintf('\n') 'Perl says: ' raw_signal];
    ensure = [ensure sprintf('\n') 'Your sequence was: ' S];
    error(ensure);
end