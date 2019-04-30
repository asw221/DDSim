
## KPD Deceased Donor Simulations

Files can be downloaded and installed on a Unix-based platform with
`make`. To get this to work you may or may not have to edit your
`~/.bash_profile`. I added the following line to mine on my Mac:

```Shell
export GUROBI_HOME='/Library/gurobi800/mac64'
```


If you have access to the group Box, data needed to run the simulation
are available [here](https://umich.app.box.com/folder/75240512306
"MBox link"). These should be extracted into a `./data` directory.

The basic command to run the simulations follows an `./DDSim.exe
path/to/parameterFile` idiom, e.g.

```Shell
usr$ ./DDSim.exe Test/Test.txt
```

The program assumes the `parameterFile` path will lie within the
`./parameters` directory; the contents of `parameterFile` should
follow a `#lowercase-varname=value` idiom, e.g.

```
## Possible contents of Test.txt
#initkpdsize=240
#maxcyclesize=3
#maxchainlength=3
#probpairinactivetoactive=0.02
#probinitkpdattrition=0.02
#probpairattrition=0.005
#probbridgedonorattrition=0.05
```



All of the simulation parameters that are settable via a
`parameterFile` and their default values are as follows:

```C++
  //Output Folder
  outputFolder = "Test";
  subFolder = "Test";

  //Simulation Settings
  utilityScheme = UTILITY_TRANSPLANTS;
	
  //Numerical Parameters
  numberOfIterations= 200;
  startingIterationID = 1;

  initKPDSize = 200;
  initKPDPairsAddedPerMatchRun = 30;
  pairArrivals = 500.0;
  nddArrivals = 10.0;

  timeSimulation = 1095; // 3 Years = 365 * 3 Days
  timeBetweenMatchRuns = 1;
  timeBetweenSelectionAndTransplantation = 30;
	
  maxCycleSize = 3;
  maxChainLength = 3;

  matchUtilityLowerBound = 1.0;
  matchUtilityUpperBound = 1.0;

  probPairActiveToInactive = 0.01;
  probPairInactiveToActive = 0.02;

  probInitKPDAttrition = 0.02;
  probPairAttrition = 0.005;
  probBridgeDonorAttrition = 0.05;
  probNodeFailure = 0.1;
  probEdgeFailure = 0.1;

  praEligibilityMin = 58;
  praEligibilityMax = 98;
	
  //Additional Options
  estimateExpectedUtility = false;
  numberOfExpectedUtilityIterations = 100;

  reserveODonorsForOCandidates = false;
  allowABBridgeDonors = false;
  allowDesensitization = true;
	
  //Files and Folders
  fileKPDData = "APDData.csv";
  fileHLAFrequency = "HLAFrequency.csv";
  fileHLADictionary = "HLADictionary.csv";
  fileSurvivalParameters = "SurvivalParameters.csv";
  fileDeceasedDonors = "DeceasedDonors.csv";
  fileWaitingListCandidates = "CandidateWaitlist.csv";
	
  //Random Number Generators Seeds
  rngSeedSelection = 9007900;
  rngSeedAttrition = 52531;
  rngSeedArrival = 5416162;
  rngSeedMatch = 3102156;
  rngSeedDonor = 3942252;
  rngSeedStatus = 7156372;
  rngSeedExpectedUtility = 1923323;
	
  //Output Suppression
  suppressExchangeInformation = true;
  suppressSimulationInformation = true;
  suppressPopulationInformation = false;
```

