function impetigounity(folder, times, varargin)

    seth = 1; setv = 1;

    if length(varargin) > 0
        if varargin{1} == 0
            setv = 0;
        end
    elseif length(varargin) < 1
        seth = 0;
    end;
    
    if length(varargin) == 2
        if varargin{2} == 0
            setv = 0;
        end
    elseif length(varargin) < 2
        setv = 0;
    end;

    d = 0; singleton = 0;
    if isdir(folder)
        d = [dir([folder '/*.txt']); dir([folder '/*.fasta'])];
    else
        d = which(folder);
        singleton = 1;
    end
    
    subfolder = 'bargraphs';
    mkdir(fullfile(folder, subfolder));

    function plot_average(file)
    
        [folder, file, ext] = fileparts(file);

        displacement = walrus_surprise([file ext]);
        global shoals sands;

        waitress = zeros(times, 0);
        for i = 1:times
            [x, waits] = displacement({}, {});
            waitress(i, 1:length(waits)) = waits;
        end
        


        tip = mean(waitress);
        
        if seth
            maxh = varargin{1};
        else
            maxh = length(tip);
        end
        
        h = figure(1); set(h, 'Renderer', 'OpenGL');
            plot(0,0); bar(1:maxh,tip(1:maxh)); title(file);
            if setv
                axis([0 maxh+1 0 varargin{2}]);
            else
                ylimits = ylim;
                axis([0 maxh+1 ylim]);
            end
            grid; xlabel('Codon Number'); ylabel('Wait Time');
        saveas(h, fullfile(folder, subfolder, [file '_bar.png']), 'png');
    
    end
    
    if singleton
        plot_average(d);
        return;
    end
    
    for i = 1:length(d)
        name = d(i).name;
        plot_average(fullfile(folder, name));
    end
       
end