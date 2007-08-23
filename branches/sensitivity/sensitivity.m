function [yields] = sensitivity(file, limit)

yields = [];

for i = -100:30
    displacement = walrus_surprise(file);
    global shoals sands;
    for j = 1:limit
        displacement([9], i*pi/180);
    end
    yields = [yields shoals/sands];
    disp(['Yield: ' num2str(shoals/sands)]);
    disp(['phi_sp: ' num2str(i)]);
    disp(' ');
end

plot(-100:30, yields);