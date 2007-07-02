function scattershot(protein, folder, leader, limit)

for i=1:limit
    seq = pearl('dysentery.pl', ['"' protein '" "' leader '"']);
    fid = fopen([folder '/' num2str(i) '.txt'], 'w');
    fprintf(fid, seq);
    fclose(fid);
end