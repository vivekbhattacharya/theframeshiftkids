function check_god(folder);
    load TAV.mat;
    TAV = TAV(find(TAV));
    norm_TAV = TAV ./ sum(TAV);

    disp('Looping through TAV vectors');
    
    d = dir([folder '/*.mat']);
    for i = 1:length(d)
        file = fullfile(folder, d(i).name);
        
        % file contains an `optimal` variable.
        load(which(file));
        norm_opt = optimal ./ sum(optimal);
        diff = abs((norm_TAV - norm_opt) ./ norm_TAV);
        fprintf('%s : %g\n', file, mean(diff));
    end
end