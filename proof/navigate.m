function [new] = navigate(char, edit, index, prfb)
index = str2num(index);
if (char == 'd')
    set(edit, 'String', disco(prfb, index, -3));
    new = index - 3;
elseif (char == 'f')
    set(edit, 'String', disco(prfb, index, -1));
    new = index - 1;
elseif (char == 'j')
    set(edit, 'String', disco(prfb, index, 1));
    new = index + 1;
elseif (char == 'k')
    set(edit, 'String', disco(prfb, index, 3));
    new = index + 3;
elseif (char == 'z')
    set(edit, 'String', disco(prfb, index, -index + 1));
    new = 1;
else
    new = index;
end

function [excerpt] = disco(prfb, index, delta)
excerpt = [prfb(index+delta : index+delta+40) ' {{' prfb(index+delta+41:index+delta+43) '}} ' prfb(index+delta+44:index+delta+90)];