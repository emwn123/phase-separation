# Default case

Reproduce the data and figures in `default/` with:

```sh
phase-separation --output reproduced-default
```

The original run used Windows, Python 3.14 and NumPy 2.5.4. Parameters are saved in `default/config.json`. Mean composition drift at saved times was below 8e-16. Free energy decreased from 1033.3885622605567 to 274.703312237177, with no increases between saved samples. This is a numerical check for this case, not a convergence result.
