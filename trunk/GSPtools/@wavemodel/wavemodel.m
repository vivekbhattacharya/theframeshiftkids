function self = wavemodel(file, fs)
    config;
    global Config;

    self  = struct();
    super = model(file, fs);
    self  = class(self, 'wavemodel', super);
end
