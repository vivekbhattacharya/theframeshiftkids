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

         h = figure(1); set(h, 'Visible', 'off');
         mesh(x, y, yields');
         title(file); grid;
         xlabel('Initial Displacement');
         ylabel('Species Angle (Degrees)');
         zlabel('Error-Free Rate');
         saveas(h, image, 'png');
         directory = fileparts(image); 
         save(fullfile(directory, [file '.mat']), 'x', 'y', 'yields');
     end
end

function [x y yields] = grope(d, file, limit)
    init_disps = -0.5:0.1:1.5;
    angles = -180:10:90;
    yields = zeros(length(init_disps), length(angles));

    global shoals sands Config;

    for i = 1:length(init_disps)
        Config.init_disp = i;
        init_disp = init_disps(i);
        
        for k = 1:length(angles)
            Config.phi_sp = k*pi/180;
            angle = angles(k);

            shoals = 0; sands = 0;
            for j = 1:limit, d([25]); end
            yields(i, j) = shoals/sands;

            disp(['  Yield: ' num2str(shoals/sands)]);
            disp(['  phi_sp: ' num2str(angle)]);
            disp(['  init_disp: ' num2str(init_disp)]);
            disp(' ');
        end
    end
end