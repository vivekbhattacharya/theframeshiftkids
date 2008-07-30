function [Travel] = load_travel()
    config; global Config;
    load(superwhich(Config.Travel));
end
