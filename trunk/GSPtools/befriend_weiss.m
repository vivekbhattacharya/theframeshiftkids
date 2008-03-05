% Look for the limit on modifying Travel.uua until we see a difference
% between 4p2003 (which should not frameshift properly but does anyway
% right now) and prfB (which is our golden boy).
function befriend_weiss(folder, limit)
     multipliers = [9.7:0.01:9.9];
     classify(folder, 'befriend_weiss', @helper);

     function helper(displacement, n, file, image)
         disp(file);
         [a, b] = fileparts(image);
         prefix = fullfile(a, b);

         [yields] = grope(multipliers, displacement, file, limit);
         disp([multipliers; yields]);

         h = figure(1); % set(h, 'Visible', 'off');
         plot(multipliers, yields);
         title(file);
         xlabel('Multipliers'); ylabel('EFR');
         saveas(h, image, 'png');
     end
end

function [yields] = grope(multipliers, d, file, limit)
    yields = zeros(1, length(multipliers));
    global shoals sands Config;

    % Abort displacement fast.
    for k = 1:length(multipliers)
        global Travel;
        Travel.uua = Travel.cuu * multipliers(k);
        shoals = 0; sands = 0;
        for j = 1:limit, d([25]); end
        yields(k) = shoals/sands;
    end
end