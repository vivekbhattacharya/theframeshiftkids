function [Travel] = load_travel()
    global Config; config;
    load(which(Config.TAV));
      Travel.aga = 500;
%     allcodons = fieldnames(Travel);
%     for i = 1:length(allcodons),
%         eval(['Travel.' allcodons{i} ' = 2 * Travel.' allcodons{i} ';']);
%     end;
%     Travel.uua = Travel.cuu * 4.7;
end