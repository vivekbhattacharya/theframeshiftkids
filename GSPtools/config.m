function config
    global Config;
    Config = struct();

    % 0: Deviation. 1: Probability.
    Config.yield = 1;

    Config.init_disp = 0.1;
    Config.phi_sp = 3 * pi/180;
    Config.signal_shift = [];

    % Normally, don't abort displacement at the first sign of trouble.
    Config.dire = 0;
    Config.should_cache = 1;

    Config.detect_pauses = 0;
    
    % Set the temperature (Celsius) and the values to use ('Freier' or 'XiaMathews')
    Config.temp = 37;
    Config.values = 'Freier';
end
