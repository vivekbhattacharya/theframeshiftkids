function sensitivity(folder, limit)
     classify(folder, 'sensitivity', @helper);
         
     function helper(displacement, n, file, image)
         disp(file);
         [x y yields] = grope(displacement, file, limit);

         h = figure(1); set(h, 'Visible', 'off');
         plot3(x, y, yields);
         title(file);
         grid;
         xlabel('Initial displacement');
         ylabel('Phase angle');
         zlabel('Error-free rate');
         saveas(h, image, 'png');
     end
end

function [x y yields] = grope(d, file, limit)
    yields = [];
    x = -0.5:0.05:1.5;
    y = -180:5:90;
    myi = 0;

    global shoals sands Config;
    for i = -0.5:0.05:1.5
        myj = 0;
        myi = myi + 1;
        for k = -180:5:90
            myj = myj + 1;
            shoals = 0; sands = 0;
            Config.phi_sp = k*pi/180;
            Config.init_disp = i;
            for j = 1:limit, d([25]); end

            yields(myi,myj) = shoals/sands;
            disp(['  Yield: ' num2str(shoals/sands)]);
            disp(['  phi_sp: ' num2str(k)]);
            disp(['  c2: ' num2str(i)]);
            disp(' ');
        end
    end
    %plot(-0.2:0.05:1.5, yields);
end