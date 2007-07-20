function [seq, count, dvec] = walrus_surprise(file)
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

[signal, seq] = get_signal(file);
[mag, theta, count] = cumm_mag_phase(signal);
[dvec, theta] = diff_vectors(mag, theta, count);

figure(2); title('Cumulative phase');
    plot(0,0); polar(theta, mag);
    xlabel('Codon'); ylabel('Phase angle (deg)');