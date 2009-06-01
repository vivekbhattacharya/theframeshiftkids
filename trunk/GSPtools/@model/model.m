% Constructor for model. Takes a filename and a possible single
% location of frameshift, usually 25 for some reason.
function obj = model(file, fs)
    config;
    global Config;

    [signal, seq] = get_signal(file);
    % Free energy vector.
    obj.signal = signal;
    % mRNA sequence.
    obj.seq    = seq(13:end);
    % Place of frameshift, if there is one.
    obj.fs     = fs;
    % Number of codons, theoretically.
    obj.upper  = floor(length(obj.signal)/3);
    % tRNA values.
    obj.travel = load_travel();

    % Ribosomal shift so far.
    obj.shift  = 0;
    % Strings of +1 frameshifts.
    obj.wal    = {};
    % Strings of +1 frameshifts.
    obj.rus    = {};
    % Locations of +1 frameshifts.
    obj.nut    = [];

    % For polynomial fitting.
    obj.power   = 4;
    obj.psignal = [zeros(1, 6) obj.signal zeros(1, 6)];

    obj = class(obj, 'model');
end
