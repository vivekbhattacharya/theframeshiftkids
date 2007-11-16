function config
    global Config;
    Config = struct();
    
    % 0: Deviation. 1: Probability.
    Config.yield = 1;
    
    Config.phi_sp = -150 * pi/180;
    Config.init_disp = 0.1;
    Config.signal_shift = [0 0 0 0];
    
    % Normally, don't abort displacement at the first sign of trouble.
    Config.dire = 1;
    
    Config.detect_stops = 0;
    Config.detect_pauses = 0;
end