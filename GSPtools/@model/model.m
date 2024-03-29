% Given an mRNA file and a possible frameshift site (usually 25 for
% some reason), returns a @model instance. @model uses Dr. Bitzer's
% fourth(fifth?)-degree on a moving window of free-energy values. See
% also: @regmodel, @wavemodel.
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
