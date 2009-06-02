% Given an mRNA file and a possible frameshift site (usually 25 for
% some reason), returns a @regmodel instance. @regmodel uses a
% register-based sine model on a moving window of free-energy values
% as detailed in Ponnala et al. See also: @wavemodel, @model.
function self = regmodel(file, fs)
    config;
    global Config;

    super        = model(file, fs);
    signal       = get(super, 'signal');

    [mag, phase] = cumm_energy(signal);
    self.dvec    = inst_energy(mag, phase);
    self         = class(self, 'regmodel', super);
end
