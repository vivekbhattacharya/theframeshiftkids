% I take a work folder and a sample size as arguments and spit out 3D
% sensitivity plots into a "sensitivity" subdirectory in the work
% folder along the lines of jejunity and opportunity.
%
% sensitivity('c:\weiss', 20)
function sensitivity(folder, limit)
     classify(folder, 'sensitivity', @helper);
         
     function helper(displacement, n, file, image)
         disp(file);
         [x y yields] = grope(displacement, file, limit);

         proprietary = fullfile(fileparts(image), file);
         save([proprietary '.mat'], 'x', 'y', 'yields');
         
         h = figure(1); %set(h, 'Visible', 'off');
         mesh(x, y, yields');
         title(file);
         xlabel('Initial Displacement');
         ylabel('Species Angle (Degrees)');
         zlabel('Error-Free Rate');
         axis([xlim ylim 0 2]);
         saveas(h, image, 'png');
         saveas(h, [proprietary '.fig'], 'fig');
     end
end

function [init_disps angles yields] = grope(d, file, limit)
    init_disps = [0:0.001:0.01];
    angles = [-125:10:50];
    %    init_disps = [-0.02:0.0025:0.005];
    %    angles = [50:10:230];
    yields = zeros(length(init_disps), length(angles));

    global shoals sands Config;

    % Abort displacement fast.
    Config.dire = 1;
    for i = 1:length(init_disps)
        Config.init_disp = init_disps(i);
        
        for k = 1:length(angles)
            angle = angles(k);
            Config.phi_sp = angle*pi/180;

            shoals = 0; sands = 0;
            for j = 1:limit, d([25]); end
            yields(i, k) = shoals/sands;

            disp(['  Yield: ' num2str(shoals/sands)]);
            disp(['  phi_sp: ' num2str(angle)]);
            disp(['  init_disp: ' num2str(Config.init_disp)]);
            disp(' ');
        end
    end
end