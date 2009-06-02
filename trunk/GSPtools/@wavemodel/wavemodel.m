% Given an mRNA file and a possible frameshift site (usually 25 for
% some reason), returns a @wavemodel instance. @wavemodel uses a
% non-register model on a moving window of free-energy values. See
% also: @regmodel, @model.
function self = wavemodel(file, fs)
    config;
    global Config;

    self  = struct();
    super = model(file, fs);
    self  = class(self, 'wavemodel', super);
end
