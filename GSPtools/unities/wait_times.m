% Draws a wait time bar plot for each gene in the work folder
% (txt or FASTA), saving that plot to wait_times\gene.png
% in the work folder using a sample size of N.
%
% If given, X specifies the upper limit of the horizontal
% axis (codon number) and Y specifies the upper limit
% of the vertical axis (wait time), both with the minimum
% at zero.
%
% wait_times('c:\work folder', N, X, Y)
function wait_times(folder, times, varargin)
    l = length(varargin);
    x_bound = 0; y_bound = 0;
    if l >= 1, x_bound = varargin{1}; end;
    if l >= 2, y_bound = varargin{2}; end;

    classify(folder, 'wait_times', @plot_average);
    function plot_average(displacement, n, filename, image)
        disp(filename);
        waitress = zeros(times, 0);
        for i = 1:times
            [x, waits] = displacement([]);
            waitress(i, 1:length(waits)) = waits;
        end
        tip = mean(waitress);

        h = figure(1); set(h, 'Renderer', 'OpenGL');
            if x_bound > 0, bar(1:x_bound, tip(1:x_bound));
            else bar(1:length(tip), tip);
            end
            if y_bound > 0, axis([xlim 0 y_bound]); end

            set(h, 'Visible', 'off'); title(filename);
            grid; xlabel('Codon Number'); ylabel('Wait Time');
        saveas(h, image, 'png');
    end
end
