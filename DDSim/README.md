
## KPD Deceased Donor Simulations

Files can be downloaded and installed on a Unix-based platform with `make`.

If you have access to the group Box, data needed to run the simulation
are available [here](https://umich.app.box.com/folder/75240512306
"MBox link").

The basic command to run the simulations follows an `./DDSim.exe
path/to/parameterFile` idiom, e.g.

```Shell
./DDSim.exe Test/Test.txt
```

The program assumes the `parameterFile` path will lie within the
`./parameters` directory; the contents of `parameterFile` should
follow a `#varname=value` idiom, e.g.

```
#initkpdsize=240
#maxcyclesize=3
#maxchainlength=3
#probpairinactivetoactive=0.02
#probinitkpdattrition=0.02
#probpairattrition=0.005
#probbridgedonorattrition=0.05
```
