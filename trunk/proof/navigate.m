function navigate(char, edit)
global index prfb;

delta = 0;
switch char
   case 'd'; delta = -3;
   case 'f'; delta = -1;
   case 'j'; delta = 1;
   case 'k'; delta = 3;
   case 'z'; delta = -index + 1;
   otherwise return;
end
set(edit, 'String', disco(index+delta));


function [excerpt] = disco(b)
excerpt = [prfb(b:b+40) ' {{' prfb(b+41:b+43) '}} ' prfb(b+44:b+90)];