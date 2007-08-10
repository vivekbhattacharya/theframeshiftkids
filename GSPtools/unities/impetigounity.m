function impetigounity(file, times)

displacement = walrus_surprise(file);
global shoals sands;

waitress = zeros(times, 0);
for i = 1:times
    [x, waits] = displacement({}, {});
    waitress(i, 1:length(waits)) = waits;
end

tip = mean(waitress);
h = figure(1); set(h, 'Renderer', 'OpenGL');
    plot(0,0); plot(1:length(tip), tip);
    grid; xlabel('Codon Number'); ylabel('Wait Time');