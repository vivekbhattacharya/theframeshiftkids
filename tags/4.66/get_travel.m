% Combines Codons and a TAV matrix.
function travel = get_travel(codons, TAV)
    limit = length(codons);
    str = 'struct(';
    
    max_hund = max(TAV);
    min_hund = min(TAV(find(TAV)));
    for i = 1:limit
        k = nloopcalc(i, TAV, max_hund, min_hund);
        str = [str sprintf('''%s'', %g,', codons{i}, k)];
    end
    
    % Delete trailing comma.
    str(end) = [];
    travel = eval([str ');']);
end

function N = nloopcalc(index, TAV, max, min)
    % Stop codons should have high wait times. Manually set values for
    % E.coli since the TAV will be zero for these codons
    if TAV(index) == 0
        N = 1000;
    else
        % Compare x with Names and extract abundance ratio
        N = max/min - floor(TAV(index)/min);
    end
end