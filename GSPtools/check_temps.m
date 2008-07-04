function [yields] = check_temps(folder, limit);
    temps = [5:0.5:55];
    classify(folder, 'sensitivity', @helper);

    function helper(displacement, n, file, image)
        disp(file);
        [a, b] = fileparts(image);
        prefix = fullfile(a, b);

        [yields] = grope(temps, displacement, file, limit);

        h = figure(1); % set(h, 'Visible', 'off');
        plot(temps, yields);
        title(file);
        xlabel('Temperature (Centigrade)');
        ylabel('Yield');
        axis([0 max(temps) + 5 0 1]);
     end
end


function [yields] = grope(temps, d, file, limit)
    yields = zeros(1, length(temps));
    global shoals sands Config;

    % Abort displacement fast.
    Config.dire = 1;
    for i = 1:length(temps)
        Config.temp = temps(i);
        disp(['  Temperature: ' num2str(temps(i))]);

        shoals = 0; sands = 0;
        for j = 1:limit, d([25]); end
        yields(i) = shoals/sands;
    end
end