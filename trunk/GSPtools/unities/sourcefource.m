function sourcefource(folder)
    classify(folder, 'sourcefource!', @helper, 'preparation');
    function helper(path, filename, image)
        disp(filename);
        
        quasiforce(filename);
        h = figure(1);
        [folder, file, ext] = fileparts(path);
        saveas(h, image, 'png');
        saveas(h, [folder file '.fig'], 'fig');
    end
end