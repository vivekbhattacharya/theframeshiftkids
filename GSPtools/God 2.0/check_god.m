function check_god(folder);
    load TAV.mat;
    norm_TAV = TAV./sum(TAV);
    norm_TAV2 = (norm_TAV == 0);
    norm_TAV2 = norm_TAV + norm_TAV2;
    disp('Looping through TAV vectors.');
    
    d = [dir([folder '/*.mat'])];
    for i = 1:length(d)
        file = fullfile(folder, d(i).name);
        load(which(file));
        norm_opt = optimal./(sum(optimal));
        diffvec = abs((norm_TAV - norm_opt)./norm_TAV2);
        disp([': ' file ' : ' num2str(mean(diffvec))]);
    end
        