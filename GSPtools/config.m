function config
    global Config;
    Config = struct();

    Config.Travel = 'Travel2.mat';

    % 0: Deviation. 1: Probability.
    Config.yield = 0;

    Config.init_disp = 0.1;
    Config.phi_sp = -30 * pi/180;
    Config.signal_shift = 0;

    % Normally, don't abort displacement at the first sign of trouble.
    Config.dire = 0;
    Config.should_cache = 1;

    % Parameters for displacement.m: C1 is the proportionality constant
    % between force and displacement. Power is the exponent of the
    % sinusoidal probability model at the heart of the stochastic part
    % of the model.
    Config.c1 = 0.005;
    Config.power = 10;

    % Set the temperature (in Celsius) and the values to use
    % (Kidnap::Freier or Kidnap::XiaMathews). Freier is
    % temperature-independent because it only works for
    % room-temperature, cf. the paper.
    Config.temp = 37;
    Config.energies = 'Kidnap::Freier';
    % auuccuccacuag for E coli. More species are listed in
    % scan_brightly.pl's help text.
    Config.tail = 'auuccuccacuag';
end
