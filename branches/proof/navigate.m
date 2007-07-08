function [delta] = navigate(char, edit)
global index; delta = 0;
switch char
   case 'd'; delta = -3;
   case 'f'; delta = -1;
   case 'j'; delta = 1;
   case 'k'; delta = 3;
   case 'z'; delta = -index + 1;
   otherwise return;
end
set(edit, 'String', disco(index+delta));
index = index + delta;

function [excerpt] = disco(b)
global prfb;
excerpt = [prfb(b:b+30) ' {{' prfb(b+31:b+33) '}} ' prfb(b+36:b+45)];