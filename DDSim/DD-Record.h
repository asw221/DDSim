/* ---------------------------------------------
   DD-Record.h
   Stores random portion of specific simulation iteration
   Ensures replicability and integrity of results
   Elements are copied to simulation at beginning of each iteration
   ---------------------------------------------- */

#ifndef RECORD_H
#define RECORD_H

#include "DD-Enums-Functions.h"
#include "DD-Candidate.h"
#include "DD-Donor.h"
#include "DD-Data.h"
#include "DD-Parameters.h"
#include "DD-Node.h"
#include "DD-Match.h"
#include "DD-MatchRun.h"
#include "DD-RNG.h"

#include <cmath>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


class KPDRecord {

private:
  KPDData * kpdData;
  KPDParameters * kpdParameters;

  int currentIteration;
  int nodeIDAssignment;  // stores the integer ID of the last KPD node assigned

  int nDonorsTotal;
  int nPairsTotal;
  int nNDDsTotal;

  // Pool Information	
  std::vector<KPDNode *> kpdNodes;	
  std::vector<std::deque<KPDStatus> > kpdNodeStateTransitionMatrix;
  std::vector<std::deque<int> > kpdNodeStateTransitionTimeMatrix;
	
  std::map<int, std::map<int, std::vector<KPDMatch *> > > kpdMatches;

  std::vector<std::vector<bool> > kpdAdjacencyMatrix;
  std::vector<std::vector<bool> > kpdAdjacencyMatrixReduced;

  // Random Number Generators
  RNG rngSelection;
  RNG rngAttrition;
  RNG rngArrival;
  RNG rngMatch;
  RNG rngDonor;
  RNG rngStatus;

  // Helper Functions
  double computeMatchSuccessProb(
    const double &assumedSuccessProb, const bool &nddDonor = false
  ) const;
  
  void clearRecord();
  void generateInitialKPD();
  // ^^ For pre-simulation: time = -T -> 0
  void assembleKPD(std::vector<int> matchRunTimes);
  // ^^ For entire simulation: time = 0 -> ...
  void assignStateTransitions();
  void assignMatchProperties();
	
  std::stringstream kpdRecordLog;

  // Add output files for initial (pre-simulation) KPD pool construction
  static int initKPDIterationCount;
  std::string fileKPDOutputPop;
  std::string fileKPDOutputTransplants;
  std::ofstream outputKPDPopulationCsv;   // Everyone who ever enters the initial KPD
  std::ofstream outputKPDTransplantsCsv;  // Transplants performed


public:
  KPDRecord(KPDData * data, KPDParameters * params);
  ~KPDRecord();
	
  KPDDonor * generateDonor();
  void generateDonors(
    std::vector<KPDDonor *> & donors,
    int numberOfDonors,
    KPDCandidate * candidate
  );
  
  KPDMatch * generateMatch(
    KPDCandidate * candidate,
    KPDDonor * donor,
    KPDCrossmatch virtualCrossmatchResult,
    bool waitlist,
    bool donorIsNDD
  );

  void generateSimulationData(int iteration, std::vector<int> matchRunTimes);

  // Cloning Functions
  std::vector<KPDNode *> getNodes();
  std::vector<KPDNodeType> getNodeTypes();
  std::vector<std::deque<KPDStatus> > getKPDNodeStateTransitionMatrix();
  std::vector<std::deque<int> > getKPDNodeStateTransitionTimeMatrix();

  std::map<int, std::map<int, std::vector<KPDMatch *> > > getKPDMatches();

  std::vector<std::vector<bool> > getAdjacencyMatrix();
  std::vector<std::vector<bool> > getAdjacencyMatrixReduced();
	
  std::string getPopulationList();

  void printLog();	
};




KPDRecord::KPDRecord(KPDData * data, KPDParameters * params){
  kpdData = data;
  kpdParameters = params;

  // #################################################################
  fileKPDOutputPop = "output/" + params->getOutputFolder() + "/" +
    params->getSubFolder() + "/InitialKPDPopulation.csv";
  fileKPDOutputTransplants = "output/" + params->getOutputFolder() + "/" +
    params->getSubFolder() + "/InitialKPDTransplants.csv";
}

int KPDRecord::initKPDIterationCount = 0;


KPDRecord::~KPDRecord(){
  printLog();
}




double KPDRecord::computeMatchSuccessProb(
  const double &assumedSuccessProb, const bool &nddDonor
) const {
  const double numberOfPotentialNodeFailures = (nddDonor) ? 1.0 : 2.0;
  // const double failureProb = (1 - assumedSuccessProb) +
  //   kpdParameters->getProbEdgeFailure() +
  //   kpdParameters->getProbNodeFailure() * numberOfPotentialNodeFailures;
  // return 1.0 - ((failureProb > 1) ? 1.0 : failureProb);
  const double pEdge = 1.0 - kpdParameters->getProbEdgeFailure();
  const double pNode = 1.0 - kpdParameters->getProbNodeFailure();
  const double successProb = assumedSuccessProb * pEdge *
    std::pow(pNode, numberOfPotentialNodeFailures);
  return successProb;
};



void KPDRecord::clearRecord() {
  nodeIDAssignment = 1;
  nDonorsTotal = 0;
  nPairsTotal = 0;
  nNDDsTotal = 0;
  
  kpdNodes.clear();
  kpdNodeStateTransitionMatrix.clear();
  kpdNodeStateTransitionTimeMatrix.clear();
  kpdMatches.clear();
  kpdAdjacencyMatrix.clear();
  kpdAdjacencyMatrixReduced.clear();
}









void KPDRecord::generateInitialKPD() {
  // ###################################################################
  initKPDIterationCount++;
  if (initKPDIterationCount == 1) {
    outputKPDPopulationCsv.open(fileKPDOutputPop);
    outputKPDTransplantsCsv.open(fileKPDOutputTransplants);
    
    std::string recipientHeader = "RecipBT,RecipPRA,RecipAge,RecipSex,"
      "RecipRace,RecipDiabetes,RecipHeight,RecipWeight,RecipBMI,"
      "RecipPreviousTransplant,RecipTOD,RecipHepC,RecipInsurance";
    std::string donorHeader = "DonorBT,DonorHLA,DonorRelation,DonorAge,"
      "DonorSex,DonorRace,DonorHeight,DonorWeight,DonorBMI,"
      "DonorCigarette";
    outputKPDTransplantsCsv << "Iteration,MatchRun,DonorNodeID,NDD,"
			 << donorHeader << ","
			 << "RecipNodeID,"
			 << recipientHeader << ","
			 << "Arrangement" << "\n";
    outputKPDPopulationCsv << "Iteration,MatchRun,RecipNodeID,NodeType,"
			<< recipientHeader << "\n";
  }
  else {
    outputKPDPopulationCsv.open(fileKPDOutputPop, std::ios_base::app);
    outputKPDTransplantsCsv.open(fileKPDOutputTransplants, std::ios_base::app);
  }
  // ###################################################################
  
  kpdRecordLog << "Initializing a Mature KPD:" << std::endl;
  std::cout << "Initializing a Mature KPD..." << std::endl;
  std::cout << "  (of size ~" << kpdParameters->getInitKPDSize() << ")\n";

  int initKPDSize = kpdParameters->getInitKPDSize();
  std::vector<KPDNode *> initialNodes;
  int initialMatchRunCount = 0;
		
  while ((int) initialNodes.size() < initKPDSize) {
    initialMatchRunCount++;
    int nPairs = 0;
    KPDDonor * newNDD = generateDonor();
    KPDNode * newNDDNode = new KPDNode(nodeIDAssignment, 0, newNDD);
    initialNodes.push_back(newNDDNode);
    nodeIDAssignment++;

    // Add ~30 Pairs (depending on getInitKPDPairsAddedPerMatchRun())
    while (nPairs < kpdParameters->getInitKPDPairsAddedPerMatchRun()) {	
      std::pair<KPDCandidate *, int> candidateInfo = kpdData->drawCandidate(
        rngSelection.runif()); // Draw candidate with replacement

      KPDCandidate * candidate = candidateInfo.first;
      int numberOfDonors = candidateInfo.second;

      std::vector<KPDDonor *>  donors;
      generateDonors(donors, numberOfDonors, candidate); // Generate donors for candidate

      if (donors.size() == numberOfDonors) {
	KPDNode * newInitialNode = new KPDNode(nodeIDAssignment, 0, donors, candidate);
	initialNodes.push_back(newInitialNode);
	nodeIDAssignment++;
				
	nPairs++;

	// ###################################################################
	outputKPDPopulationCsv << currentIteration << ","
			    << initialMatchRunCount << ","
			    << newInitialNode->getID() << ","
			    << KPDFunctions::nodeTypeToString(newInitialNode->getType()) << ","
			    << newInitialNode->getCandidateString()
			    << "\n";
	// ###################################################################
      }
    }

    std::map<int, std::map<int, std::vector<KPDMatch *> > > initialMatches;
    int N = (int) initialNodes.size();
		
    for (int i = 1; i <= N; i++) {
      int donorNodeIndex = i - 1;
      KPDNode * donorNode = initialNodes[donorNodeIndex];
      int donorNodeID = donorNode->getID();

      //Get donor nodes
      std::vector<KPDDonor *> donors = donorNode->getDonors();
      int numDonors = donorNode->getNumberOfDonors();


      for (int j = 1; j <= N; j++) {
	int candidateNodeIndex = j - 1;

	KPDNode * candidateNode = initialNodes[candidateNodeIndex];
	int candidateNodeID = candidateNode->getID();

	if (i != j && candidateNode->getType() == PAIR) {
	  bool allowableMatchExists = false;
	  std::vector<KPDMatch *> matches;
	  //Iterate through associated donors
	  
	  for (int k = 1; k <= numDonors; k++) {
	    int donorIndex = k - 1;
	    //Perform virtual crossmatch
	    KPDCrossmatch virtualCrossmatchResult = kpdData->performCrossmatch(
              candidateNode->getCandidate(), donors[donorIndex], false);

	    if (kpdData->allowableMatch(virtualCrossmatchResult) &&
		(rngMatch.runif() <
		 computeMatchSuccessProb(1.0, donorNode->isNDD()))) {
	      KPDMatch * match = new KPDMatch();
	      match->setAdjacency(true);
	      match->setSuccessfulMatch(true);
	      match->setVirtualCrossmatchResult(virtualCrossmatchResult);
	      allowableMatchExists = true;
	      matches.push_back(match);
	    }
	    else {
	      KPDMatch * noMatch = new KPDMatch();
	      noMatch->setVirtualCrossmatchResult(virtualCrossmatchResult);
	      matches.push_back(noMatch);
	    }
	  }

	  if (allowableMatchExists) {
	    initialMatches[donorNodeID][candidateNodeID] = matches;
	  }
	}
      }
    }
		
    std::vector<KPDStatus> initialStatuses(N, STATUS_ACTIVE);
    std::vector<KPDTransplant> initialTransplanted(N, TRANSPLANT_NO);
    std::vector<KPDNodeType> initialNodeTypes;

    for (std::vector<KPDNode *>::iterator it = initialNodes.begin();
	 it != initialNodes.end(); it++) {
      initialNodeTypes.push_back((*it)->getType());
    }

    KPDMatchRun * initialMatchRun = new KPDMatchRun(
      kpdParameters, 0, 0, initialNodes, initialNodeTypes,
      initialStatuses, initialTransplanted, initialMatches);

    std::vector<std::vector<int> > cyclesAndChains;
    std::vector<double> utilities;

    // Find all the LRSs in the current pool
    initialMatchRun->collectCyclesAndChainsForCurrentMatchRun(cyclesAndChains);

    // If there are indeed arrangements in the pool...
    if (cyclesAndChains.size() > 0) {
      // Assign the appropriate utility values to the LRSs
      initialMatchRun->assignUtilitiesForCurrentMatchRun(cyclesAndChains, utilities);
      
      // Select optimal set of LRSs
      std::vector<int> optimalSolution;
      initialMatchRun->getOptimalSolutionForCurrentMatchRun(
        optimalSolution, cyclesAndChains, utilities);

      // std::cout << "  -- GUROBI CALL  --\n";

      for (std::vector<int>::iterator itSolution = optimalSolution.begin();
	   itSolution != optimalSolution.end(); itSolution++) {

	std::vector<int> cycleOrChain = cyclesAndChains[*itSolution];

	// #####################################################################
	// Add output to outputTransplantsCsv
	// (code could be a bit fragile)
	//
	
	bool isChain = false;
	bool donorFound = false, recipientFound = false, donorSubIndexFound = false;
	int currentRecipientID, currentDonorID, donorSubIndex;
	KPDNode* currentDonor;
	KPDNode* currentRecipient;
	for (std::vector<int>::iterator itCycle = cycleOrChain.begin();
	     itCycle != cycleOrChain.end(); itCycle++) {
	  
	  // Check current arrangement is a cycle or chain & find donor
	  for (std::vector<KPDNode *>::iterator itNode = initialNodes.begin();
	       itNode != initialNodes.end(); itNode++) {
	    if ((*itNode)->getID() == *itCycle) {
	      if (itCycle == cycleOrChain.begin()) {
		isChain = (*itNode)->isNDD();  // NDD's always come first
	      }
	      currentDonor = (*itNode);
	      currentDonorID = currentDonor->getID();
	      donorFound = true;
	    }
	  }

	  // Find recipient
	  if (donorFound) {
	    if ((itCycle - cycleOrChain.begin()) < (cycleOrChain.size() - 1)) {
	      currentRecipientID = (*(itCycle + 1));
	    }
	    else if (!isChain) {
	      currentRecipientID = cycleOrChain[0];
	    }
	    else {
	      break;  // <- a little sloppy
	    }
	    for (std::vector<KPDNode *>::iterator itNode = initialNodes.begin();
		 itNode != initialNodes.end(); itNode++) {
	      if ((*itNode)->getID() == currentRecipientID) {
		currentRecipient = (*itNode);
		recipientFound = true;
	      }
	    }

	    if (recipientFound) {  // print output
	      // std::map<int, std::map<int, std::vector<KPDMatch *> > > initialMatches;
	      // initialMatches[donorNodeID][candidateNodeID] = matches;
	      for (int donorIndex = 0;
		   donorIndex < initialMatches[currentDonorID][currentRecipientID].size();
		   donorIndex++) {
		if (initialMatches[currentDonorID][currentRecipientID][donorIndex]
		    ->getAdjacency()) {
		  donorSubIndex = donorIndex;
		  donorSubIndexFound = true;
		}
	      }

	      if (donorSubIndexFound) {
		outputKPDTransplantsCsv << currentIteration << ","
					<< initialMatchRunCount << ","
					<< currentDonorID << ","
					<< KPDFunctions::boolToTF(currentDonor->isNDD()) << ","
					<< currentDonor->getDonorString(donorSubIndex) << ","
					<< currentRecipientID << ","
					<< currentRecipient->getCandidateString()
					<< ",";
		for (std::vector<int>::iterator itc = cycleOrChain.begin();
		     itc != cycleOrChain.end(); itc++) {
		  if (itc == cycleOrChain.begin())
		    outputKPDTransplantsCsv << (*itc);
		  else
		    outputKPDTransplantsCsv << "-" << (*itc);
		}
		outputKPDTransplantsCsv << "\n";
	      }
	    }
	  }  // if (donorFound)
	}  // for (std::vector<int>::iterator itCycle = cycleOrChain.begin(); ...
	// #####################################################################
	
	for (std::vector<int>::iterator itCycle = cycleOrChain.begin();
	     itCycle != cycleOrChain.end(); itCycle++) {
	  for (std::vector<KPDNode *>::iterator itNode = initialNodes.begin();
	       itNode != initialNodes.end(); itNode++) {
	    if ((*itNode)->getID() == *itCycle) {
	      initialNodes.erase(itNode);
	      break;
	    }
	  }
	}
      }	 // for (std::vector<int>::iterator itSolution = optimalSolution.begin();
    }  // if (cyclesAndChains.size() > 0)

    
    int i = 0;
    while (i < (int) initialNodes.size()) {		
      double u = rngAttrition.runif();
      if (u < kpdParameters->getProbInitKPDAttrition()) {
	initialNodes.erase(initialNodes.begin() + i);
      }
      else {
	i++;
      }
    }
  }  // while ((int) initialNodes.size() < initKPDSize)

  
  kpdRecordLog << "Initial KPD Pool: " << std::endl;

  for (std::vector<KPDNode *>::iterator itNode = initialNodes.begin();
       itNode != initialNodes.end(); itNode++) {

    KPDNode * node = *itNode;
    kpdNodes.push_back(node);
    if (node->getType() == NDD) {
      nNDDsTotal++;
      kpdRecordLog << "NDD " << node->getID()
		   << ": " << node->getDonor(0)->donorOutput() << std::endl;
    }
    else {
      int numberOfDonors = node->getNumberOfDonors();
      nPairsTotal++;
      nDonorsTotal += numberOfDonors;			

      kpdRecordLog << "Pair " << node->getID() << std::endl;
      kpdRecordLog << node->getCandidate()->candidateOutput() << std::endl;
      for (int donorIndex = 0; donorIndex < numberOfDonors; donorIndex++) {
	kpdRecordLog << node->getDonor(donorIndex)->donorOutput() << std::endl;
      }
      kpdRecordLog << std::endl;
    }
  }

  kpdRecordLog << "Initial Pool NDDs: " << nNDDsTotal
	       << ", Pairs: " << nPairsTotal
	       << ", Donors: " << nDonorsTotal << std::endl;
  std::cout << "Initial Pool NDDs: " << nNDDsTotal
	    << ", Pairs: " << nPairsTotal
	    << ", Donors: " << nDonorsTotal << std::endl;

  // ###################################################################
  outputKPDTransplantsCsv.close();
  outputKPDPopulationCsv.close();
  // ###################################################################
}
// void KPDRecord::generateInitialKPD()







void KPDRecord::assembleKPD(std::vector<int> matchRunTimes) {
  kpdRecordLog << "Selecting Nodes:" << std::endl;
  std::cout << "Selecting Nodes..." << std::endl;
	
  double pairArrivals = kpdParameters->getPairArrivals() / 365.0;
  double nddArrivals = kpdParameters->getNDDArrivals() / 365.0;
	
  // Generate NDDs	
  double nddTimeTracker = rngArrival.rexp(nddArrivals); // Poisson Process for NDD arrival
  double pairTimeTracker = rngArrival.rexp(pairArrivals); // Poisson Process for pair arrival

  int time = 1;

  for (std::vector<int>::iterator it = matchRunTimes.begin();
       it != matchRunTimes.end(); it++) {
    int matchRunTime = *it;
    int nNDDs = 0;
    int nPairs = 0;

    while (time <= matchRunTime) {
      while (nddTimeTracker <= time) { // Generate NDDs
	KPDDonor * newNDD = generateDonor(); // Randomly generate NDD
	KPDNode * newNDDNode = new KPDNode(nodeIDAssignment, time, newNDD);
	kpdNodes.push_back(newNDDNode);
	nodeIDAssignment++;

	kpdRecordLog << "At Time " << time << ", NDD "
		     << newNDDNode->getID() << " Joins" << std::endl;
	kpdRecordLog << newNDD->donorOutput() << std::endl;
	kpdRecordLog << std::endl;
	std::cout << "At Time " << time
		  << ", NDD " << newNDDNode->getID()
		  << " Joins" << std::endl;

	nddTimeTracker += rngArrival.rexp(nddArrivals); // Get next NDD arrival time
	nNDDs++;
      }			

      while (pairTimeTracker <= time) { // Generate pairs	
	bool successfulPairGeneration = false;
	while (!successfulPairGeneration) {
	  std::pair<KPDCandidate *, int> candidateInfo =
	    kpdData->drawCandidate(rngSelection.runif()); // Draw candidate with replacement

	  KPDCandidate * candidate = candidateInfo.first;
	  int numberOfDonors = candidateInfo.second;
	  std::vector<KPDDonor *>  donors;
	  generateDonors(donors, numberOfDonors, candidate); // Generate donors for candidate

	  if ((int)donors.size() == numberOfDonors) {
	    nDonorsTotal += numberOfDonors;
	    KPDNode * newPairNode = new KPDNode(nodeIDAssignment, time, donors, candidate);
	    kpdNodes.push_back(newPairNode);
	    nodeIDAssignment++;

	    kpdRecordLog << "At Time " << time << ", Pair "
			 << newPairNode->getID() << " Joins" << std::endl;
	    kpdRecordLog << candidate->candidateOutput() << std::endl;
	    
	    for (int donorIndex = 0; donorIndex < (int)donors.size(); donorIndex++) {
	      kpdRecordLog << donors[donorIndex]->donorOutput() << std::endl;
	    }
	    kpdRecordLog << std::endl;
	    std::cout << "At Time " << time << ", Pair "
		      << newPairNode->getID() << " Joins" << std::endl;

	    pairTimeTracker += rngArrival.rexp(pairArrivals); // Get next pair arrival time
	    nPairs++;
	    successfulPairGeneration = true;
	  }
	}
      }
      time++;
    }  // while (time <= matchRunTime)

    kpdRecordLog << "At Match Run at Time " << matchRunTime
		 << ", # of NDDs Added: " << nNDDs
		 << ", # of Pairs Added: " << nPairs
		 << std::endl;
    std::cout << "At Match Run at Time " << matchRunTime
	      << ", # of NDDs Added: " << nNDDs
	      << ", # of Pairs Added: " << nPairs
	      << std::endl;
    nNDDsTotal += nNDDs;
    nPairsTotal += nPairs;
  }  // for (std::vector<int>::iterator it = matchRunTimes.begin(); ...

  kpdRecordLog << "Total Pairs: " << nPairsTotal
	       <<  ", Total NDDs: " << nNDDsTotal
	       << ", Total Donors: " << nDonorsTotal
	       << std::endl;
  std::cout << "Total Pairs: " << nPairsTotal
	    << ", Total NDDs: " << nNDDsTotal
	    << ", Total Donors: " << nDonorsTotal
	    << std::endl;

}
// void KPDRecord::assembleKPD(std::vector<int> matchRunTimes)








void KPDRecord::assignStateTransitions() {	
  kpdRecordLog << "Generating State Transitions:" << std::endl;
  std::cout << "Generating State Transitions..." << std::endl;
	
  int timeSimulation = kpdParameters->getTimeSimulation();
  int timeBetweenSelectionAndTransplantation =
    kpdParameters->getTimeBetweenSelectionAndTransplantation();

  double probPairAttrition = kpdParameters->getProbPairAttrition();
  double probBridgeDonorAttrition = kpdParameters->getProbBridgeDonorAttrition();
  double probPairActiveToInactive = kpdParameters->getProbPairActiveToInactive();
  double probPairInactiveToActive = kpdParameters->getProbPairInactiveToActive();
		
  //Iterate through KPD nodes
  for (std::vector<KPDNode *>::iterator it = kpdNodes.begin();
       it != kpdNodes.end(); it++) {

    KPDNode * node = *it;
    int nodeID = node->getID();

    kpdRecordLog << "Node " << nodeID << std::endl;

    std::deque<KPDStatus> kpdNodeStateTransitions;
    std::deque<int> kpdNodeStateTransitionTimes;

    int time = node->getArrivalTime();
    KPDStatus status = STATUS_ACTIVE;

    kpdNodeStateTransitions.push_back(status);
    kpdNodeStateTransitionTimes.push_back(time);

    //State transitions for pairs
    if (node->getType() == PAIR) {
      
      kpdRecordLog << nodeID << ": " << time
		   << " (" << KPDFunctions::statusToString(status) << ") ";
					
      while (time < (timeSimulation + timeBetweenSelectionAndTransplantation) &&
	     status != STATUS_WITHDRAWN) {			
	time++;
	double w = rngStatus.runif(); // Withdrawal
	double u = rngStatus.runif(); // Active to Inactive or Inactive to Active
	if (w < probPairAttrition) {
	  status = STATUS_WITHDRAWN;
	  kpdNodeStateTransitions.push_back(status);
	  kpdNodeStateTransitionTimes.push_back(time);
	  kpdRecordLog << time << " (" << KPDFunctions::statusToString(status) << ") ";
	}
	else if (status == STATUS_ACTIVE && u < probPairActiveToInactive) {
	  status = STATUS_INACTIVE;
	  kpdNodeStateTransitions.push_back(status);
	  kpdNodeStateTransitionTimes.push_back(time);
	  kpdRecordLog << time << " (" << KPDFunctions::statusToString(status) << ") ";
	}
	else if (status == STATUS_INACTIVE && u < probPairInactiveToActive) {
	  status = STATUS_ACTIVE;
	  kpdNodeStateTransitions.push_back(status);
	  kpdNodeStateTransitionTimes.push_back(time);
	  kpdRecordLog << time << " (" << KPDFunctions::statusToString(status) << ") ";
	}
      }
      kpdRecordLog << std::endl;
    }
    else {
      int nddDeparts = time + 90;
      
      kpdRecordLog << nodeID << ": " << time
		   << " (" << KPDFunctions::statusToString(status) << ") ";
      
      kpdNodeStateTransitions.push_back(STATUS_WITHDRAWN);
      kpdNodeStateTransitionTimes.push_back(nddDeparts);

      kpdRecordLog << nddDeparts
		   << " (" << KPDFunctions::statusToString(STATUS_WITHDRAWN)
		   << ") " << std::endl;
    }

    kpdNodeStateTransitionMatrix.push_back(kpdNodeStateTransitions);
    kpdNodeStateTransitionTimeMatrix.push_back(kpdNodeStateTransitionTimes);
  }
}








void KPDRecord::assignMatchProperties() {
  kpdRecordLog << "Generating Match Properties:" << std::endl;
  std::cout << "Generating Match Properties..." << std::endl;

  int N = (int) kpdNodes.size();

  kpdAdjacencyMatrix.assign(1 + N, std::vector<bool>(1 + N, false));
  kpdAdjacencyMatrixReduced.assign(1 + N, std::vector<bool>(1 + N, false));

  //Iterate through donor nodes
  for (int i = 1; i <= N; i++) {
    int donorNodeIndex = i - 1;

    KPDNode * donorNode = kpdNodes[donorNodeIndex];
    int donorNodeID = donorNode->getID();

    //Get donor nodes
    std::vector<KPDDonor *> donors = donorNode->getDonors();
    int numDonors = donorNode->getNumberOfDonors();

    //Iterate through candidate nodes
    for (int j = 1; j <= N; j++) {
      int candidateNodeIndex = j - 1;
      KPDNode * candidateNode = kpdNodes[candidateNodeIndex];
      int candidateNodeID = candidateNode->getID();
			
      if (i != j) {
	// Pair -> NDD (Implicit Backward Edge from all Donors to the NDD)
	if ((donorNode->getType() == PAIR) &&
	    (candidateNode->getType() == NDD))
	  kpdAdjacencyMatrix[i][j] = true;

	// Pair
	else if (candidateNode->getType() == PAIR) {
	  double pra = candidateNode->getCandidatePRA();					
	  bool allowableMatchExists = false;
	  std::vector<KPDMatch *> matches;

	  //Iterate through associated donors
	  for (int k = 1; k <= numDonors; k++) {
	    int donorIndex = k - 1;

	    //Perform virtual crossmatch
	    KPDCrossmatch virtualCrossmatchResult = kpdData
	      ->performCrossmatch(candidateNode->getCandidate(),
				  donors[donorIndex], false);
					
	    if (kpdData->allowableMatch(virtualCrossmatchResult)) {
	      allowableMatchExists = true;
	      kpdAdjacencyMatrix[i][j] = true;
	      kpdAdjacencyMatrixReduced[i][j] = true;

	      KPDMatch * newMatch = generateMatch(candidateNode->getCandidate(),
						  donorNode->getDonor(donorIndex),
						  virtualCrossmatchResult, false,
						  donorNode->isNDD());
							
	      matches.push_back(newMatch);
							
	      kpdRecordLog << "   " << donorNodeID << "[" << k
			   << "] -> " << candidateNodeID << " "
			   << newMatch->matchShortOutput() << std::endl;
	    }

	    else {

	      KPDMatch * noMatch = new KPDMatch();
	      noMatch->setVirtualCrossmatchResult(virtualCrossmatchResult);

	      matches.push_back(noMatch);
	    }
	  }

	  //std::cout << donorNodeID << "->" << candidateNodeID << std::endl;
	  if (allowableMatchExists)
	    kpdMatches[donorNodeID][candidateNodeID] = matches;
	}
      }
    }
  }
}









KPDDonor * KPDRecord::generateDonor() {
  std::vector<double> u;
  for (int i = 1; i <= 5; i++) {
    u.push_back(rngDonor.runif());
  }
  KPDDonor * newDonor = kpdData->generateDonor(u);
  return newDonor;
}





void KPDRecord::generateDonors(
  std::vector<KPDDonor *> & donors,
  int numberOfDonors,
  KPDCandidate * candidate
) {

  int generatedDonors = 0;
  int attempts = 0;
	
  while (attempts < 2 * numberOfDonors && generatedDonors < numberOfDonors) {
    // This could take a while...
    KPDDonor * donor = generateDonor();		
    KPDCrossmatch crossmatch = kpdData->performCrossmatch(candidate, donor, false);

    if (!kpdData->allowableMatch(crossmatch)) {
      donors.push_back(donor);
      generatedDonors++;
    }
    attempts++;
  }
}






KPDMatch * KPDRecord::generateMatch(
  KPDCandidate * candidate,
  KPDDonor * donor,
  KPDCrossmatch virtualCrossmatchResult,
  bool waitlist,
  bool donorIsNDD
) {
  //Assign utility values
  double fiveYearSurvival = kpdData->calculateSurvival(candidate, donor, 1);
  double tenYearSurvival = kpdData->calculateSurvival(candidate, donor, 0);

  double kdpi = donor->getKDPI();
  int pra = candidate->getPRA();
  double epts = candidate->getEPTS();
  // std::cout << "EPTS: " << epts << "\n";

  double transplantDifficultyScore = 0.0001;
  if (pra >= 97 || donor->getBT() == BT_AB) {
    transplantDifficultyScore = 1.0;
  }

  double assignedUtility = rngMatch.runif(
    kpdParameters->getMatchUtilityLowerBound(),
    kpdParameters->getMatchUtilityUpperBound());

  //Generate probability values

  double assumedMatchSuccessProbability = 1.0;
  double actualMatchSuccessProbability = 1.0;

  // Set success rates before accounting for additional friction on edges/nodes
  if (waitlist) {
    // Values from TABLE in writeup. For EPTS raw -> percentile translations
    // see OPTN document in 'man' directory
    //
    if (epts < 1.51597) {                     // epts <= 20%
      if (kdpi < 0.2) {
	assumedMatchSuccessProbability = 0.9;
	actualMatchSuccessProbability = 0.9;
      }
      if (0.2 <= kdpi && kdpi < 0.5) {
	assumedMatchSuccessProbability = 0.6;
	actualMatchSuccessProbability = 0.6;
      }
      if (0.5 <= kdpi && kdpi < 0.75) {
	assumedMatchSuccessProbability = 0.3;
	actualMatchSuccessProbability = 0.3;
      }
      if (0.75 <= kdpi && kdpi < 1.0) {
	assumedMatchSuccessProbability = 0.3;
	actualMatchSuccessProbability = 0.3;
      }
    }  // if (epts < 1.51597)

    if (1.51597 <= epts && epts < 2.21624) {  // epts \in [21% - 50%]
      if (kdpi < 0.2) {
	assumedMatchSuccessProbability = 0.9;
	actualMatchSuccessProbability = 0.9;
      }
      if (0.2 <= kdpi && kdpi < 0.5) {
	assumedMatchSuccessProbability = 0.9;
	actualMatchSuccessProbability = 0.9;
      }
      if (0.5 <= kdpi && kdpi < 0.75) {
	assumedMatchSuccessProbability = 0.9;
	actualMatchSuccessProbability = 0.9;
      }
      if (0.75 <= kdpi && kdpi < 1.0) {
	assumedMatchSuccessProbability = 0.9;
	actualMatchSuccessProbability = 0.9;
      }
    }  // if (1.51597 <= epts && epts < 2.21624)

    if (2.21625 <= epts && epts < 2.61330) {  // epts \in [51% - 75%]
      if (kdpi < 0.2) {
	assumedMatchSuccessProbability = 0.0;
	actualMatchSuccessProbability = 0.0;
      }
      if (0.2 <= kdpi && kdpi < 0.5) {
	assumedMatchSuccessProbability = 0.5;
	actualMatchSuccessProbability = 0.5;
      }
      if (0.5 <= kdpi && kdpi < 0.75) {
	assumedMatchSuccessProbability = 0.5;
	actualMatchSuccessProbability = 0.5;
      }
      if (0.75 <= kdpi && kdpi < 1.0) {
	assumedMatchSuccessProbability = 1.0;
	actualMatchSuccessProbability = 1.0;
      }
    }  // if (2.0 <= epts && epts < 3.0)

    if (2.61330 <= epts) {                    // epts > 75%
      if (kdpi < 0.2) {
	assumedMatchSuccessProbability = 0.0;
	actualMatchSuccessProbability = 0.0;
      }
      if (0.2 <= kdpi && kdpi < 0.5) {
	assumedMatchSuccessProbability = 0.0;
	actualMatchSuccessProbability = 0.0;
      }
      if (0.5 <= kdpi && kdpi < 0.75) {
	assumedMatchSuccessProbability = 0.1;
	actualMatchSuccessProbability = 0.1;
      }
      if (0.75 <= kdpi && kdpi < 1.0) {
	assumedMatchSuccessProbability = 0.2;
	actualMatchSuccessProbability = 0.2;
      }
    }  // if (3.0 <= epts)
    
  }  
  else {		// i.e. else ==> !waitlist
    if (pra < 25) {
      assumedMatchSuccessProbability = 0.95;
      actualMatchSuccessProbability = 0.95;
    }
    else if (pra >= 25 && pra < 50) {
      assumedMatchSuccessProbability = 0.8;
      actualMatchSuccessProbability = 0.8;
    }
    else if (pra >= 50 && pra < 75) {
      assumedMatchSuccessProbability = 0.65;
      actualMatchSuccessProbability = 0.65;
    }
    else {
      assumedMatchSuccessProbability = 0.5;
      actualMatchSuccessProbability = 0.5;
    }
  }

  // bool successfulMatch = rngMatch.runif() < actualMatchSuccessProbability;
  // Include additional friction on edges/nodes
  const bool successfulMatch =
    rngMatch.runif() < computeMatchSuccessProb(
      assumedMatchSuccessProbability, donorIsNDD);

  KPDMatch * newMatch = new KPDMatch(
    true, fiveYearSurvival, tenYearSurvival,
    transplantDifficultyScore, assignedUtility, 
    assumedMatchSuccessProbability, actualMatchSuccessProbability,
    virtualCrossmatchResult, successfulMatch);
	
  return newMatch;
}







void KPDRecord::generateSimulationData(int iteration, std::vector<int> matchRunTimes) {
  currentIteration = iteration;
	
  kpdRecordLog << "--------------" << std::endl;
  kpdRecordLog << "Iteration " << iteration << std::endl;
  kpdRecordLog << "--------------" << std::endl << std::endl;
	
  //Initialize random number seeds

  clearRecord();

  kpdRecordLog << "Initializing Random Number Generator Seeds" << std::endl;
  std::cout << "Initializing Random Number Generator Seeds... " << std::endl;

  rngSelection.setSeed(kpdParameters->getRNGSeedSelection() * iteration);
  rngAttrition.setSeed(kpdParameters->getRNGSeedAttrition() * iteration);
  rngArrival.setSeed(kpdParameters->getRNGSeedArrival() * iteration);
  rngMatch.setSeed(kpdParameters->getRNGSeedMatch() * iteration);
  rngDonor.setSeed(kpdParameters->getRNGSeedDonor() * iteration);
  rngStatus.setSeed(kpdParameters->getRNGSeedStatus() * iteration);

  //Select new simulation data, assign node and match properties
  generateInitialKPD();
  assembleKPD(matchRunTimes);
  assignStateTransitions();
  assignMatchProperties();

  kpdRecordLog << "Simulation Data for Iteration " << iteration
	       << " Generated!" << std::endl << std::endl;
  std::cout << "Simulation Data for Iteration " << iteration
	    << " Generated!" << std::endl << std::endl;
}








std::vector<KPDNode *> KPDRecord::getNodes() {
  std::vector<KPDNode *> nodes; 
  for (std::vector<KPDNode*>::iterator it = kpdNodes.begin();
       it != kpdNodes.end(); it++)
    nodes.push_back((*it)->copy()); // Deep Copy
  return nodes;
}








std::vector<KPDNodeType> KPDRecord::getNodeTypes() {
  std::vector<KPDNodeType> nodeTypes;
  for (std::vector<KPDNode*>::iterator it = kpdNodes.begin(); it != kpdNodes.end(); it++) {
    nodeTypes.push_back((*it)->getType());
  }
  return nodeTypes;
}







std::vector<std::vector<bool> > KPDRecord::getAdjacencyMatrix() {
  std::vector<std::vector<bool> > adjacencyMatrixClone(kpdAdjacencyMatrix);
  return adjacencyMatrixClone;
}





std::vector<std::vector<bool> > KPDRecord::getAdjacencyMatrixReduced() {
  std::vector<std::vector<bool> > adjacencyMatrixReducedClone(
    kpdAdjacencyMatrixReduced);
  return adjacencyMatrixReducedClone;
}





std::vector<std::deque<KPDStatus> > KPDRecord::getKPDNodeStateTransitionMatrix() {
  std::vector<std::deque<KPDStatus> > kpdNodeStateTransitionMatrixClone(
    kpdNodeStateTransitionMatrix);
  return kpdNodeStateTransitionMatrixClone;
}






std::vector<std::deque<int> > KPDRecord::getKPDNodeStateTransitionTimeMatrix() {
  std::vector<std::deque<int> > kpdNodeStateTransitionTimeMatrixClone(
    kpdNodeStateTransitionTimeMatrix);
  return kpdNodeStateTransitionTimeMatrixClone;
}






std::map<int, std::map<int, std::vector<KPDMatch*> > > KPDRecord::getKPDMatches() {
  std::map<int, std::map<int, std::vector<KPDMatch*> > > kpdMatchesCopy;
  
  for (std::map<int, std::map<int, std::vector<KPDMatch*> > >::iterator
	 itDonor = kpdMatches.begin(); itDonor != kpdMatches.end(); itDonor++) {
    int donorNodeID = itDonor->first;
    std::map<int, std::vector<KPDMatch*> > donorMatches = itDonor->second;
    
    for (std::map<int, std::vector<KPDMatch *> >::iterator
	   itCandidate = donorMatches.begin();
	 itCandidate != donorMatches.end(); itCandidate++) {
      int candidateNodeID = itCandidate->first;
      std::vector<KPDMatch*> matches = itCandidate->second;
      std::vector<KPDMatch *> matchCopy;
      
      for (std::vector<KPDMatch *>::iterator
	     itMatch = matches.begin(); itMatch != matches.end(); itMatch++)
	matchCopy.push_back((*itMatch)->copy());
      
      kpdMatchesCopy[donorNodeID][candidateNodeID] = matchCopy;
    }
  }
  return kpdMatchesCopy;
};






std::string KPDRecord::getPopulationList() {
  std::stringstream population;
  //Iterates through each pair involved in the simulation and collects
  // demographic information
  for (std::vector<KPDNode *>::iterator itNode = kpdNodes.begin();
       itNode != kpdNodes.end(); itNode++){

    int numberOfDonors = (*itNode)->getNumberOfDonors();
    for (int donorIndex = 0; donorIndex < numberOfDonors; donorIndex++){
      population << currentIteration << ","
		 << (*itNode)->getID() << ","
		 << (donorIndex + 1) << ","
		 << KPDFunctions::nodeTypeToString((*itNode)->getType()) << ","
		 << (*itNode)->getArrivalTime() << ","
		 << (*itNode)->getCandidateString() << ","
		 << (*itNode)->getDonorString(donorIndex) << ","
		 << "\n";
    }
  }

  return population.str();
}




void KPDRecord::printLog(){	
  std::string logFile = "output/" + kpdParameters->getOutputFolder() + "/" +
    kpdParameters->getSubFolder() + "/Log-Record.txt";
  std::ofstream outputFileLog(logFile.c_str());
  outputFileLog << kpdRecordLog.str();
  outputFileLog.close();
}


#endif
