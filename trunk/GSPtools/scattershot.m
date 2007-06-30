function scattershot(protein, folder, limit)

for i=1:limit
    seq = pearl('dysentery.pl', ['"' protein '"']);
    fid = fopen([folder '/' num2str(i) '.txt'], 'w');
    fprintf(fid, seq);
    fclose(fid);
end