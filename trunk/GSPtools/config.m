function config
    global Config;
    Config = struct();
    
    % 0: Deviation. 1: Probability.
    Config.yield = 0;
    
    Config.phi_sp = -30 * pi / 180;
end