function config
    global Config;
    Config = struct();
    
    % 0: Deviation. 1: Probability.
    Config.yield = 1;
    
    Config.phi_sp = -30 * pi / 180;
    Config.init_disp = 0.005;
    
    % Normally, don't abort displacement at the first sign of trouble.
    Config.dire = 0;
    
    Config.detect_stops = 0;
    Config.detect_pauses = 0;
end