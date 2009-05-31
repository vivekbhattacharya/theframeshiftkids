function [energy, polyforce] = polyenergy(codon)
    global store;
    new_signal = [zeros(1, 6) store.signal zeros(1, 6)];
    center = 3*codon + 5 + store.shift;

    power = 4;
    lower = floor(-power/2);
    upper = floor(power/2);
    x = lower:upper;
    y = new_signal(center+lower:center+upper);
    energy = polyfit(x, y, power);

    function [dx, force] = helper(codon, x)
        force = -energy .* (power:-1:0);
        force = force(1:4);
        dx    = polyval(force, x/2 - store.shift);
    end
    polyforce = @helper;
end

