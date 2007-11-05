% ------------------------------------------------
% Plots displacment error bars for gene.txt, saving
% the plot to jejunity\gene.png in the work folder
% using an N sample size for each gene (txt or
% FASTA) in the work folder.
% 
% Usage: jejunity('C:\work folder', N)
% ------------------------------------------------
function jejunity(folder, limit)
    classify(folder, 'jejunity', @helper);
    
    function helper(displacement, n, file, image)
        disp(file);

        upper = ceil(n);
        a = zeros(limit, upper);
        for i=1:limit
            x = displacement([]);
            maxiderm = length(x);
            for j=1:upper
                if j <= maxiderm, a(i,j) = x(j);
                else a(i,j) = -3.14;
                end;
            end
        end
        [avg, stddev] = exmean(a);
        
        max_dom = 50*(ceil(length(x)/50));
        h = figure(1); set(h, 'Renderer', 'OpenGL');
            set(h, 'Visible', 'off');
            errorbar(1:upper, avg, stddev); title(file);
            
            grid; xlabel('Codon Number'); ylabel('x(k)');
            axis([0, max_dom, min(-3, floor(min(x))), max(3, ceil(max(x)))]);
        saveas(h, image, 'png');
    end
end