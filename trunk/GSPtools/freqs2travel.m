% ARGUMENTS
%
% freqs_file: Path to a slightly modified codon frequencies file with
%   a codon and its frequency separated by a space on each line.
%
% Names: From the original, the sexy, the one-and-only Codons.mat. I
%   use it to order the final TAV.
function [travel] = freqs2tav(freqs_file, Names)
    freqs_file = fopen(freqs_file);
    results = textscan(freqs_file, '%s %f');
    fclose(freqs_file);

    codons = lower(results{1});
    freqs = num2cell(results{2}/sum(results{2}));
    freq = cell2struct(freqs', codons', 2);
    tav = cell2struct(freqs', codons', 2);

    add = {
        {'gcu' 'gcu' 'gcc'},
        {'gcc' 'gcu' 'gcc'},
        {'gca' 'gca' 'gcg'},
        {'gcg' 'gca' 'gcg'},
        {'aga' 'aga' 'agg'},
        {'agg' 'aga' 'agg'},
        {'cga' 'cga' 'cgg'},
        {'cgg' 'cga' 'cgg'},
        {'cgu' 'cgu' 'cgc'},
        {'cgc' 'cgu' 'cgc'},
        {'gau' 'gau' 'gac'},
        {'gac' 'gau' 'gac'},
        {'aau' 'aau' 'aac'},
        {'aac' 'aau' 'aac'},
        {'ugu' 'ugu' 'ugc'},
        {'ugc' 'ugu' 'ugc'},
        {'gaa' 'gaa' 'gag'},
        {'gag' 'gaa' 'gag'},
        {'caa' 'caa' 'cag'},
        {'cag' 'caa' 'cag'},
        {'gga' 'gga' 'ggg'},
        {'ggg' 'gga' 'ggg'},
        {'ggu' 'ggu' 'ggc'},
        {'ggc' 'ggu' 'ggc'},
        {'cau' 'cau' 'cac'},
        {'cac' 'cau' 'cac'},
        {'auu' 'auu' 'auc'},
        {'auc' 'auu' 'auc'},
        {'uua' 'uua' 'uug'},
        {'uug' 'uua' 'uug'},
        {'cua' 'cua' 'cug'},
        {'cug' 'cua' 'cug'},
        {'cuu' 'cuu' 'cuc'},
        {'cuc' 'cuu' 'cuc'},
        {'aaa' 'aaa' 'aag'},
        {'aag' 'aaa' 'aag'},
        {'uuu' 'uuu' 'uuc'},
        {'uuc' 'uuu' 'uuc'},
        {'cca' 'cca' 'ccg'},
        {'ccg' 'cca' 'ccg'},
        {'ccu' 'ccu' 'ccc'},
        {'ccc' 'ccu' 'ccc'},
        {'agu' 'agu' 'agc'},
        {'agc' 'agu' 'agc'},
        {'uca' 'uca' 'ucg'},
        {'ucg' 'uca' 'ucg'},
        {'ucu' 'ucu' 'ucc'},
        {'ucc' 'ucu' 'ucc'},
        {'aca' 'aca' 'acg'},
        {'acg' 'aca' 'acg'},
        {'uau' 'uau' 'uac'},
        {'uac' 'uau' 'uac'},
        {'gua' 'gua' 'gug'},
        {'gug' 'gua' 'gug'},
        {'guu' 'guu' 'guc'},
        {'guc' 'guu' 'guc'},
          };

    for i = 1:length(add)
        line = add{i};
        tav.(line{1}) = freq.(line{2}) + freq.(line{3});
    end

    % Stop codons set to 0. get_travel will set it to 1000.
    tav.uga = 0; tav.uag = 0; tav.uaa = 0;

    final_tav = zeros(1, length(tav));
    for i = 1:length(Names)
        final_tav(i) = tav.(Names{i});
    end

    travel = get_travel(Names, final_tav);
end
