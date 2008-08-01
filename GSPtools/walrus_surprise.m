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
function [fantastic, colon_spastic, n] = walrus_surprise(file, varargin)
    clear global sands shoals Config;
    config();

    [signal, seq] = get_signal(file);
    aforce = force(signal);
    n = length(aforce) - 1;

    function [x, waits] = helper(fs)
        [x, waits] = displacement(seq(13:end), aforce, fs);
    end
    fantastic = @helper;

    % Must be figure(2) for alopeciaunity to work
    function draw_polar(fs)
    % Empty pending discussion.
    end
    colon_spastic = @draw_polar;

    clear globals;
    global Travel Names;
    Travel = load_travel();
end
