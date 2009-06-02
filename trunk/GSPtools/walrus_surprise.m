% So you're walking to your car at Sea World when
% you hear the chirps and groans of a large sea
% mammal flying at 900 knots toward your head and
% then and only then do you realize that you have
% been nailed by a
%
%       WALRUS SURPRISE.
%       Walrus Surprise.
%       walrus surprise.
%       WALRUS SURPRISE.
%
% Arguments: A file of genes
% Returns: Function that yields displacement,
% function that creates a polar plot, and the
% number of nucleotides.
function [fantastic, n] = walrus_surprise(file, varargin)
    clear global sands shoals Config;
    config();

    [signal, seq] = get_signal(file);
    n = floor(length(signal)/3);

    function [x] = helper(fs)
        x = displacement(seq(13:end), signal, fs);
    end
    fantastic = @helper;

    clear global;
    global Travel Names;
    Travel = load_travel();
end
