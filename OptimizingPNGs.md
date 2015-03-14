Say you've called `opportunity('C:\proteins', 10)` and now you have a folder `C:\proteins\thepictureshow\` brimming with PNG files. Did you know you can halve the folder space?

  * Download [OptiPNG](http://optipng.sourceforge.net/) and place that in your `C:\Windows`.
  * Open Command Prompt and `cd C:\proteins\thepictureshow`.
  * `for %f in (*.png) do optipng "%f"`