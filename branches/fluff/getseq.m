% GETSEQ
% This function extracts the sequence from the sequence-file
% Usage: S = getseq(SEQFILE)
% 
% SEQFILE = file containing only the sequence portion spread out over many
% lines. Each line starts with the index of the first base. For example,
% 
% Contents of SEQFILE:
%         1 agcttttcat tctgactgca acgggcaata tgtctctgtg tggattaaaa aaagagtgtc
%        61 tgatagcagc ttctgaactg gttacctgcc gtgagtaaat taaaatttta ttgacttagg
%       121 tcactaaata ctttaaccaa tataggcata gcgcacagac agataaaaat tacagagtac
%       181 acaacatcca tgaaacgcat tagcaccacc attaccacca ccatcaccat taccacaggt
%       241 aacggtgcgg gctgacgcgt acaggaaaca cagaaaaaag cccgcacctg acagtgcggg
%       301 cttttttttt cgaccaaagg taacgaggta acaaccatgc gagtgttgaa gttcggcggt
%       361 acatcagtgg caaatgcaga acgttttctg cgtgttgccg atattctgga aagcaatgcc
%       421 aggcaggggc aggtggccac cgtcctctct gcccccgcca aaatcaccaa ccacctggtg
%       .....................................................................
%       .....................................................................
%       4639621 gcatgatatt gaaaaaaata tcaccaaata aaaaacgcct tagta
% 
% WARNING: SEQFILE must be a full path. Use `which` to find full paths of
% relative paths.

function S = getseq(SEQFILE)
[status, S] = system(['perl "' which('getseq.pl') '" "' SEQFILE '"']);
