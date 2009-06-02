function val = get(self, name)
    switch name
      case 'shift'
        val = self.shift;
      case 'psignal'
        val = self.psignal;
      otherwise
        error([name ' is not a valid @model property'])
    end
end
