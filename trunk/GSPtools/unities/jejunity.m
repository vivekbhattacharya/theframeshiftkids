% ------------------------------------------------
% Plots displacment error bars for gene.txt, saving
% the plot to carnivale\gene.png in the work folder
% using an N sample size for each gene (txt or
% FASTA) in the work folder.
% 
% Usage: jejunity('C:\work folder', N)
% ------------------------------------------------
function jejunity(folder, limit)
    classify(folder, 'carnivale', @helper);
    
    function helper(displacement, n, file, image)
        disp(file);

        upper = ceil(n);
        a = zeros(limit, upper);
        for i=1:limit
            x = displacement({}, {});
            maxiderm = length(x);
            for j=1:upper
                if j <= maxiderm, a(i,j) = x(j);
                else a(i,j) = -3.14;
                end;
            end
        end
        [avg, stddev] = exmean(a);
        
        top = 3; bottom = -3;
        if max(x) > 3, top = ceil(max(x)); end;
        if min(x) < -3, bottom = floor(min(x)); end;
        
        max_dom = 50*(ceil(length(x)/50));
        h = figure(1); set(h, 'Renderer', 'OpenGL');
            set(h, 'Visible', 'off');
            plot(0,0); errorbar(1:upper, avg, stddev); title(file);
            grid; xlabel('Codon Number'); ylabel('x(k)');
            axis([0, max_dom, bottom, top]);
        saveas(h, image, 'png');
    end
end