% I take a work folder and a sample size as arguments and spit out 3D
% sensitivity plots into a "sensitivity" subdirectory in the work
% folder along the lines of jejunity and opportunity.
%
% sensitivity('c:\weiss', 20)
function sensitivity(folder, limit)
     init_disps = [-1:0.2:1.8];
     angles = [-100:10:100];
     classify(folder, 'sensitivity', @helper);
         
     function helper(displacement, n, file, image)
         disp(file);
         [yields] = grope(init_disps, angles, displacement, file, limit);

         proprietary = fullfile(fileparts(image), file);
         save([proprietary '.mat'], 'init_disps', 'angles', 'yields');
         
         h = figure(1); set(h, 'Visible', 'off');
         mesh(init_disps, angles, yields');
         title(file);
         xlabel('Initial Displacement');
         ylabel('Species Angle (Degrees)');
         zlabel('Error-Free Rate');
         axis([xlim ylim 0 2]);
         saveas(h, image, 'png');
         saveas(h, [proprietary '.fig'], 'fig');
     end
end

function [yields] = grope(init_disps, angles, d, file, limit)
    yields = zeros(length(init_disps), length(angles));
    global shoals sands Config;

    % Abort displacement fast.
    Config.dire = 1;
    for i = 1:length(init_disps)
        Config.init_disp = init_disps(i);
        disp(['  init_disp: ' num2str(init_disps(i))]);
        
        for k = 1:length(angles)
            angle = angles(k);
            Config.phi_sp = angle*pi/180;

            shoals = 0; sands = 0;
            for j = 1:limit, d([25]); end
            yields(i, k) = shoals/sands;
        end
        disp([angles; yields(i, 1:length(angles))]);
        disp ' ';
    end
end