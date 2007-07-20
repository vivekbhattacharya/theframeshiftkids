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
% Returns: Those genes, how many there were, and
%   a vector of differences
function fantastic = walrus_surprise(file)
    [signal, seq] = get_signal(file);
    [mag, theta] = cumm_mag_phase(signal);
    [dvec, theta] = diff_vectors(mag, theta);
    
    function [x] = helper(fs, bs)
        x = displacement(seq(13:end), dvec, fs, bs);
    end
    fantastic = @helper;

    return;
    figure(2); title('Cumulative phase');
        plot(0,0); polar(theta, mag);
        xlabel('Codon'); ylabel('Phase angle (deg)');
end