# FORBILD Head phantom
The source code is for MATLAB, and it  was taken from the paper:

[Yu, Z., Noo, F., Dennerlein, F., Wunderlich, A., Lauritsch, G., & 
Hornegger, J. (2012). Simulation tools for two-dimensional experiments in 
x-ray computed tomography using the FORBILD head phantom. 
Physics in medicine and biology, 57(13),
N237.](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3426508/)

## How to use it
To generate a simple phantom of resolution 512x512, one should write:
```
phantom = forbild_gen(512, 512, 1, 1)
```

This phantom will include both left and right ears, which is
specified by flags 1 and 1 respectively.

To save generated phantom as a TIFF file, one may use:
```
convertDataTo32Tiff('path-to-dir/filename.tif', phantom)
```


