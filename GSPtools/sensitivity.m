% I take a work folder and a sample size as arguments and spit out 3D
% sensitivity plots into a "sensitivity" subdirectory in the work
% folder along the lines of jejunity and opportunity.
%
% sensitivity('c:\weiss', 20)
function sensitivity(folder, limit)
     powers = [18:2:24];
     c1s = [0.000425:.000005:0.00049];
     classify(folder, 'sensitivity', @helper);

     function helper(displacement, n, file, image)
         disp(file);
         [a, b] = fileparts(image);
         prefix = fullfile(a, b);

         [yields] = grope(powers, c1s, displacement, file, limit);
         save([prefix '.mat'], 'powers', 'c1s', 'yields');

%          h = figure(1); % set(h, 'Visible', 'off');
%          mesh(powers, c1s, yields');
%          title(file);
%          xlabel('Initial Displacement');
%          ylabel('Species Angle (Degrees)');
%          zlabel('Error-Free Rate');
%          axis([xlim ylim 0 2]);
%          saveas(h, image, 'png');
%          saveas(h, [prefix '.fig'], 'fig');
     end
end

function [yields] = grope(powers, c1s, d, file, limit)
    yields = zeros(length(powers), length(c1s));
    global shoals sands Config;

    % Abort displacement fast.
    Config.dire = 1;
    for i = 1:length(powers)
        Config.power = powers(i);
        disp(['  Power: ' num2str(powers(i))]);

        for k = 1:length(c1s)
            Config.c1 = c1s(k)*pi/180;

            shoals = 0; sands = 0;
            for j = 1:limit, d([25]); end
            yields(i, k) = shoals/sands;
        end
        disp([c1s; yields(i, 1:length(c1s))]);
        disp ' ';
    end
end