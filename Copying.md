The first step is to copy the SVN repository to a different place,
which will preserve the ancient and beautiful history of our code.
You'll need an intermediate directory on your local hard drive to do
so. For example, "C:/repos/BWFtools-copy". To do so, open the command
prompt and type

```
> cd C:/repos/
> svnadmin create BWFtools-copy
> svnsync init file:///repos/BWFtools-copy http://theframeshiftkids.googlecode.com/svn/trunk/GSPtools
> synsync sync file:///repos/BWFtools-copy
[long time]
```

Now suppose you opened a Google Code Project Hosting account and your
HTTPS SVN URL was "https://me.googlecode.com/svn/trunk". (You could
stop there, but you could also stick this code in a subdirectory, such
as "http://me.googlecode.com/svn/trunk/my-tools".) A note about
naming: BWFtools is an unregistered trademark, and you cannot use it
in creating a name for your project.

To copy your temporary repository over to the one on Google, do

```
> cd C:/repos/
> svnsync init -u username https://me.../trunk file:///repos/BWFtools-copy
> svnsync sync -u username https://me.../trunk
[long time]
```

where "username" is your Google Code username. The full process
[is documented by Google](http://code.google.com/support/bin/answer.py?answer=56673&topic=10386).

The second thing you'll have to do is add this license to all code
files (or in a "COPYING" file under the trunk or on a web page, even
if neither of these two options are kosher).

```
Original BWFtools code copyright (c) 2007, Hao Lian and Vivek Bhattacharya
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the
      distribution.

    * Neither the name of the <organization> nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

This software is provided by the copyright holders and contributors
"as is" and any express or implied warranties, including, but not
limited to, the implied warranties of merchantability and fitness for
a particular purpose are disclaimed. In no event shall the copyright
owner or contributors be liable for any direct, indirect, incidental,
special, exemplary, or consequential damages (including, but not
limited to, procurement of substitute goods or services; loss of use,
data, or profits; or business interruption) however caused and on any
theory of liability, whether in contract, strict liability, or tort
(including negligence or otherwise) arising in any way out of the use
of this software, even if advised of the possibility of such damage.
```

You also must retain `_COPYING` under "/pearls/Kidnap" and its original
text in some form accessible to all project contributors and users.