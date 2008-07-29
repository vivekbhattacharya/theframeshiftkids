function config
    global Config;
    Config = struct();

    Config.TAV = 'Travel2.mat';

    % 0: Deviation. 1: Probability.
    Config.yield = 0;

    Config.init_disp = 0.1;
    Config.phi_sp = -30 * pi/180;
    Config.signal_shift = [];

    % Normally, don't abort displacement at the first sign of trouble.
    Config.dire = 0;
    Config.should_cache = 1;

    Config.detect_pauses = 0;

    % Parameters for displacement.m
    Config.c1 = 0.005;
    Config.power = 10;

    % Set the temperature (Celsius) and the values to use ('Freier' or 'XiaMathews')
    Config.temp = 37;
    Config.values = 'Freier';
    Config.tail = 'auuccuccacuag'; % auuccuccacuag for E coli
end
