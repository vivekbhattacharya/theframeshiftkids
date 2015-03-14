A quick note: Many functions below take folders as as arguments. However, if you would rather run one file, you might do this: Move the file to its own folder and then pass that argument. However, you can pass just a file. For example, below, you can type `yields('prfB.txt', [25], 5)` into Matlab and it'll work the same. In fact, it'll function as if you _had_ moved the file to another folder. But you didn't.

## yields.m ##
```
M> yields('C:\up\genes\ecoli', [], 5);
aaaD.txt: 0.262499 +/- 0.0643838
aaeA.txt: 1.1049 +/- 0.148593
aaeB.txt: 1.06465 +/- 0.0372803
```

A common scenario that arises when tinkering with the model is assessing the yields for a folder of genes. Rather than running megaunity on each gene, there's `yields.m`. It takes three arguments: a folder, an array of frameshift sites, and the number of times to sample the yield. Remember, the yield is defined by `Config.yield` in `config.m`. For deviations, which have little variation, a size of 5 suffices. For probabilities, use a size of 25 instead. A larger sample size creates a smaller confidence interval but takes longer to run.