/* ---------------------------------------------
DD-Arrangement.h
Stores Found Exchanges, Define Functions to Find
Fallback Options, Calculate Utility
---------------------------------------------- */

#ifndef ARRANGEMENT_H
#define ARRANGEMENT_H

#include "DD-Enums-Functions.h"
#include "DD-Node.h"
#include "DD-Match.h"
#include "DD-Parameters.h"

#include <vector>
#include <string>
#include <deque>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "gurobi_c++.h"



class KPDArrangement {

private:
	
  //Arrangement Information
  int currentIteration;
  int _tempApproxPoolSize;
  double matchRunTime;
  double transplantationTime;
  double assignedUtility;

  bool arrangementHasNDD;

  std::vector<KPDNode *> arrangementNodes;
  std::vector<KPDNodeType> arrangementNodeTypes;

  std::vector<bool> arrangementAvailability;
	
  std::map<int, std::map<int, std::vector<KPDMatch *> > > arrangementMatches;

  std::vector<std::vector<int> > arrangementTransplantSelected;

  std::vector<std::vector<int> > arrangementFallbackOptions;
  std::vector<double> arrangementUtilityOfFallbackOptions;

  
  //Helper Functions
  
  int getChild(
    int lower,
    int current,
    std::vector<int> & visited,
    std::vector<std::vector<bool> > & adjacency
  );
  
  // indexOf() returns the integer index of nodeID's position in
  // arrangementMatches
  int indexOf(int nodeID);
	

public:
  KPDArrangement(
    int iteration,
    double mrTime,
    double txTime,
    double util,
    std::vector<KPDNode *> arrangement,
    std::vector<KPDNodeType> nodeTypes,
    std::map<int, std::map<int, std::vector<KPDMatch *> > > &matches
  );

  ~KPDArrangement();
	
  bool hasNDD();
  bool isChain();

  std::vector<int> getNodeIDs();
  int getNumberOfNodes();

  double getTransplantationTime();

  bool getAdjacency(
    int donorNodeID,
    int candidateNodeID,
    int donorIndex
  );
  
  bool getAvailability(int nodeID);

  void setFailure(int nodeID);

  int selectDonor(
    int donorNodeID,
    int candidateNodeID,
    KPDUtilityScheme utilityScheme,
    bool external
  );

  void getOptimalSolutionForFallBackOptions(
    std::vector<std::vector<int> > & options,
    int maxCycleSize,
    int maxChainLength,
    bool allowABBridgeDonors,
    KPDUtilityScheme utilityScheme
  );	

  std::string toLabel(const std::string &_delim = "-");
  std::string toTransplantString();
};




KPDArrangement::KPDArrangement(
  int iteration,
  double mrTime,
  double txTime,
  double util,
  std::vector<KPDNode *> arrangement,
  std::vector<KPDNodeType> nodeTypes,
  std::map<int, std::map<int, std::vector<KPDMatch *> > > &matches  // do not touch!
) {
  _tempApproxPoolSize = matches.size();
  int N = (int)arrangement.size();

  //Copy Arrangement
  currentIteration = iteration;
  matchRunTime = mrTime;
  transplantationTime = txTime;
  assignedUtility = util;

  arrangementHasNDD = false;

  for (int i = 1; i <= N; i++){
    int arrangementIndex = i - 1;
    arrangementNodes.push_back(arrangement[arrangementIndex]);
    arrangementNodeTypes.push_back(nodeTypes[arrangementIndex]);
    if (nodeTypes[arrangementIndex] != PAIR) {
      arrangementHasNDD = true;
    }
  }
	
  //Initialize Arrangement Properties

  arrangementAvailability.assign(N, true); // All candidates need to be available
  arrangementTransplantSelected.assign(N, std::vector<int>(N, -1));

  //Assign Arrangement Properties
  for (int i = 1; i <= N; i++) {
    int donorNodeIndex = i - 1;
    int donorNodeID = arrangementNodes[donorNodeIndex]->getID();
    int numDonors = arrangementNodes[donorNodeIndex]->getNumberOfDonors();

    for (int j = 1; j <= N; j++) {
      int candidateNodeIndex = j - 1;
      int candidateNodeID = arrangementNodes[candidateNodeIndex]->getID();

      if (matches.find(donorNodeID) != matches.end()) {
	if (matches[donorNodeID].find(candidateNodeID) != matches[donorNodeID].end()) {
	  for (int k = 1; k <= numDonors; k++) {
	    int donorIndex = k - 1;
	    arrangementMatches[donorNodeIndex][candidateNodeIndex]
	      .push_back(matches[donorNodeID][candidateNodeID][donorIndex]->copy());
	  }
	}  // if (matches[donorNodeID].find(candidateNodeID) != matches[donorNodeID].end())
      }  // if (matches.find(donorNodeID) != matches.end())
    }  // for (int j = 1; j <= N; j++)
  }  // for (int i = 1; i <= N; i++)
}
// KPDArrangement::KPDArrangement(...)




KPDArrangement::~KPDArrangement(){
  arrangementNodes.clear();
  arrangementNodeTypes.clear();
  arrangementAvailability.clear();
  arrangementMatches.clear();
  arrangementTransplantSelected.clear();
  arrangementFallbackOptions.clear();
  arrangementUtilityOfFallbackOptions.clear();
}




int KPDArrangement::getChild(
  int lower,
  int current,
  std::vector<int> & visitedVector,
  std::vector<std::vector<bool> > & adjacency
) {
  int nV = (int)visitedVector.size() - 1;
  for (int j = lower + 1; j <= nV; j++) {
    if (adjacency[current][j] == true && visitedVector[j] == 0)
      return j;
  }
  return -1;
}





int KPDArrangement::indexOf(int nodeID) {
  for (int i = 0; i < (int)arrangementNodes.size(); i++) {
    if (arrangementNodes[i]->getID() == nodeID) {
      return i;
    }
  }
  return -1;
}




bool KPDArrangement::hasNDD() {
  return arrangementHasNDD;
}


bool KPDArrangement::isChain() {
  // return arrangementHasNDD;
  bool chain = false;
  std::vector<KPDNodeType>::iterator it = arrangementNodeTypes.begin();
  while (it != arrangementNodeTypes.end() && !chain) {
    chain = *it != PAIR;
    it++;
  }
  return chain;
};




std::vector<int> KPDArrangement::getNodeIDs() {
  std::vector<int> nodeIDs;
  for (std::vector<KPDNode *>::iterator nodeIt = arrangementNodes.begin();
       nodeIt != arrangementNodes.end(); nodeIt++) {
    nodeIDs.push_back((*nodeIt)->getID());
  }
  return nodeIDs;
}




int KPDArrangement::getNumberOfNodes(){
  return (int)arrangementNodes.size();
}




double KPDArrangement::getTransplantationTime() {
  return transplantationTime;
}




bool KPDArrangement::getAdjacency(int donorNodeID, int candidateNodeID, int donorIndex) {
  bool adjacent = false;
  int donorNodeIndex = indexOf(donorNodeID);
  int candidateNodeIndex = indexOf(candidateNodeID);
  
  if (arrangementMatches.find(donorNodeIndex) != arrangementMatches.end()) {
    if (arrangementMatches[donorNodeIndex].find(candidateNodeIndex) !=
	arrangementMatches[donorNodeIndex].end()) {
      if (donorIndex < 0 ||
	  donorIndex >= arrangementMatches[donorNodeIndex][candidateNodeIndex].size()) {
	std::logic_error ("KPDArrangement::getAdjacency : donorIndex out of bounds");
      }
      adjacent = arrangementMatches[donorNodeIndex][candidateNodeIndex]
	[donorIndex]->getAdjacency();
    }
  }
  // return false;   // AW - ?????
  return adjacent;
}



bool KPDArrangement::getAvailability(int candidateID) {
  return arrangementAvailability[indexOf(candidateID)];
}



void KPDArrangement::setFailure(int nodeID){
  const int nodeIndex = indexOf(nodeID);
  if (nodeIndex > 0 && nodeIndex < arrangementAvailability.size()) {
    arrangementAvailability[nodeIndex] = false;
  }
  else {
    std::cout << "  <Check KPDArrangement::setFailure for nodeID "
	      << nodeID << ">\n";
    // This flag is usually harmless. Often this routine is called
    // without checking if the node has been transplanted/removed
    // (in which case setFailure() will ping with the message).
    // Included for completeness
  }
}






int KPDArrangement::selectDonor(
  int donorNodeID,
  int candidateNodeID,
  KPDUtilityScheme utilityScheme,
  bool external
) {
  int donorNodeIndex;
  int candidateNodeIndex;
  if (external){
    donorNodeIndex = indexOf(donorNodeID);
    candidateNodeIndex = indexOf(candidateNodeID);
  }
  else {
    donorNodeIndex = donorNodeID;
    candidateNodeIndex = candidateNodeID;
  }

  double maxVal = 0.0;
  int maxIndex = -1; // No Donor Found
  int numDonors = arrangementNodes[donorNodeIndex]->getNumberOfDonors();
  

  for (int k = 1; k <= numDonors; k++) {
    int donorIndex = k - 1;
    if (arrangementMatches.find(donorNodeIndex) != arrangementMatches.end() &&
	(arrangementMatches[donorNodeIndex].find(candidateNodeIndex) !=
	 arrangementMatches[donorNodeIndex].end())) {
      if (arrangementMatches[donorNodeIndex][candidateNodeIndex][donorIndex]
	  ->getAdjacency()) {

	if (utilityScheme == UTILITY_TRANSPLANTS) {
	  if (maxVal == 0.0) {
	    maxVal = 1.0;
	    maxIndex = donorIndex;
	  }
	}
	else if (arrangementMatches[donorNodeIndex][candidateNodeIndex]
		 [donorIndex]->getUtility(utilityScheme) > maxVal) {
	  maxVal = arrangementMatches[donorNodeIndex][candidateNodeIndex]
	    [donorIndex]->getUtility(utilityScheme);
	  maxIndex = donorIndex;
	}			
      }
    }
  }

  if (external) {
    arrangementTransplantSelected[donorNodeIndex][candidateNodeIndex] = maxIndex;
  }

  return maxIndex;
}









void KPDArrangement::getOptimalSolutionForFallBackOptions(
  std::vector<std::vector<int> > & options,
  int maxCycleSize,
  int maxChainLength,
  bool allowABBridgeDonors,
  KPDUtilityScheme utilityScheme
) {

  int maximum = std::max(maxCycleSize, maxChainLength + 1);

  // Build (small, N x N) adjacency matrix -------
  int N = (int)arrangementNodes.size();
  std::vector<std::vector<bool> > adjacencyMatrix(1 + N, std::vector<bool>(1 + N, false));

  for (int i = 1; i <= N; i++) {
    int donorNodeIndex = i - 1;
    int numDonors = arrangementNodes[donorNodeIndex]->getNumberOfDonors();

    for (int j = 1; j <= N; j++) {
      int candidateNodeIndex = j - 1;
      if (i != j) {
	if (arrangementAvailability[candidateNodeIndex]) {
	  if (arrangementNodeTypes[donorNodeIndex] == PAIR &&
	      arrangementNodeTypes[candidateNodeIndex] != PAIR) {
	    adjacencyMatrix[i][j] = true;
	  }
	  else if (arrangementNodeTypes[candidateNodeIndex] == PAIR) {
	    if (arrangementMatches.find(donorNodeIndex) !=
		arrangementMatches.end()) {
	      
	      if (arrangementMatches[donorNodeIndex].find(candidateNodeIndex) !=
		  arrangementMatches[donorNodeIndex].end()) {

		for (int k = 1; k <= numDonors; k++) {
		  int donorIndex = k - 1;
		  if (arrangementMatches[donorNodeIndex][candidateNodeIndex]
		      [donorIndex]->getAdjacency() &&
		      arrangementMatches[donorNodeIndex][candidateNodeIndex]
		      [donorIndex]->getSuccessfulMatch()) {
		    adjacencyMatrix[i][j] = true;
		    break;
		  }
		}
	      }  // if (arrangementMatches[donorNodeIndex].find(...
	    }  // if (arrangementMatches.find(donorNodeIndex) !=
	  }  // if (arrangementNodeTypes[donorNodeIndex] == PAIR && ...
	}  // if (arrangementAvailability[candidateNodeIndex])
      }  // if (i != j)
    }  // for (int j = 1; j <= N; j++)
  }  // for (int i = 1; i <= N; i++)

  // Add implicit back-edges
  for (int j = 0; j < N; j++) {
    bool isCandidateNDD = arrangementNodeTypes[j] != PAIR;
    if (isCandidateNDD) {
      for (int i = 0; i < N; i++)
	if (i != j)  adjacencyMatrix[i + 1][j + 1] = true;
    }
  }
	
  // ---------------------------------------------

  
  //Find remaining options
  int start = 1;
  int nV = (int)arrangementNodes.size();
  std::vector<int> visitedVector(nV + 1, 0);
  std::vector<int> stack_vec;

  while (start <= nV) {
    visitedVector[start] = 1;

    stack_vec.push_back(start);
    int v = getChild(start, start, visitedVector, adjacencyMatrix);
    while (!stack_vec.empty()) {
      if (v == -1) {
	int top = stack_vec.back();
	stack_vec.pop_back();
	if (top == start) {
	  visitedVector[top] = 0;
	  break;
	}
	visitedVector[top] = 0;
	v = getChild(top, stack_vec.back(), visitedVector, adjacencyMatrix);
      }
      else {
	visitedVector[v] = 1;
	stack_vec.push_back(v);
	
	if (adjacencyMatrix[v][start]) {
	  int multipleNDDCheck = 0;
	  int index = 0;
	  std::vector<int> potentialCycleOrChain;

	  for (int i = 0; i < (int)stack_vec.size(); i++) {
	    potentialCycleOrChain.push_back(stack_vec[i] - 1);
	  }

	  for (int i = 0; i < (int)potentialCycleOrChain.size(); i++) {
	    if (arrangementNodeTypes[potentialCycleOrChain[i]] != PAIR) {
	      multipleNDDCheck++;
	      index = i;
	    }
	  }

	  //Chain
	  if (multipleNDDCheck == 1 && potentialCycleOrChain.size() <= maxChainLength + 1) {
	    //Queue up NDD to front of chain
	    while (index > 0) {
	      int temp = *(potentialCycleOrChain.begin());
	      potentialCycleOrChain.erase(potentialCycleOrChain.begin());
	      potentialCycleOrChain.push_back(temp);
	      index--;
	    }

	    if (!allowABBridgeDonors) {
	      bool hasNonABBridgeDonor = false;
	      for (int k = 1; k <= arrangementNodes[*(potentialCycleOrChain.end() - 1)]
		     ->getNumberOfDonors(); k++) {

		int donorIndex = k - 1;
		if (arrangementNodes[*(potentialCycleOrChain.end() - 1)]
		    ->getDonorBT(donorIndex) != BT_AB) {
		  hasNonABBridgeDonor = true;
		  break;
		}
	      }

	      if (hasNonABBridgeDonor) {
		arrangementFallbackOptions.push_back(potentialCycleOrChain);
	      }
	    }
	    else {
	      arrangementFallbackOptions.push_back(potentialCycleOrChain);
	    }
	  }

	  //Cycle
	  else if (multipleNDDCheck == 0 && potentialCycleOrChain.size() <= maxCycleSize) {
	    arrangementFallbackOptions.push_back(potentialCycleOrChain);
	  }
	}

	if ((int)stack_vec.size() >= maximum)
	  v = -1;
	else
	  v = getChild(start, v, visitedVector, adjacencyMatrix);
      }
    }
    start++;
  }


  int numOptions = (int)arrangementFallbackOptions.size();
  if (numOptions == 0) return;
  
	
  for (std::vector<std::vector<int> >::iterator
	 itOption = arrangementFallbackOptions.begin();
       itOption != arrangementFallbackOptions.end(); itOption++) {
		
    double u = 0;

    for (std::vector<int>::iterator itNode = (*itOption).begin();
	 itNode != (*itOption).end() - 1; itNode++) {

      int donorNodeIndex = *itNode;
      int candidateNodeIndex = *(itNode + 1);
      int donorIndex = selectDonor(donorNodeIndex, candidateNodeIndex, utilityScheme, false);

      if (utilityScheme == UTILITY_TRANSPLANTS) { // Transplants

	if (arrangementNodeTypes[candidateNodeIndex] == PAIR) {	  
	  u += 1;
	}
      }
      else {	
	if (arrangementNodeTypes[candidateNodeIndex] == PAIR)
	  u += arrangementMatches[donorNodeIndex][candidateNodeIndex]
	    [donorIndex]->getUtility(utilityScheme);
      }
    }  // for (std::vector<int>::iterator itNode = (*itOption).begin() ...

    
    // Wrap from HERE in an if (!isChain) {} clause
    if (!isChain()) {
      int donorNodeIndex = *(itOption->end() - 1);
      int candidateNodeIndex = *(itOption->begin());
      int donorIndex = selectDonor(donorNodeIndex, candidateNodeIndex, utilityScheme, false);
      // (Back edge donor)

      if (utilityScheme == UTILITY_TRANSPLANTS) { // Transplants
	if (arrangementNodeTypes[candidateNodeIndex] == PAIR) {
	  u += 1;
	}
      }
      else {
	u += arrangementMatches[donorNodeIndex][candidateNodeIndex]
	  [donorIndex]->getUtility(utilityScheme);
      }
    }
    // To HERE
    
    arrangementUtilityOfFallbackOptions.push_back(u);
  }


	
  std::vector<int> selectedFallbacks;

  GRBEnv *env = 0; //Created outside try; needs to be recollected;
  GRBVar *vars = 0; //Created outside try; needs to be recollected;
  GRBVar *myVars = new GRBVar[arrangementFallbackOptions.size() + 1];
  //Created outside try; needs to be recollected;

  try{
    env = new GRBEnv();
    env->set(GRB_IntParam_Threads, 4); // Limit the threads

    GRBModel model = GRBModel(*env);
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);// 0 means no solver output; 1 otherwise;

    for (int i = 1; i <= numOptions; i++){
      std::string s = "cycle_" + KPDFunctions::intToString(i);
      myVars[i] = model.addVar(0.0, 1.0,
        -1 * arrangementUtilityOfFallbackOptions[i - 1], GRB_BINARY, s);
      //The default is min, hence the -1
    }
    model.update(); //Lazy update;

    for (int i = 0; i <= (int)arrangementNodes.size(); i++){
      GRBLinExpr expr = 0;
      for (int j = 0; j <= numOptions - 1; j++){
	for (std::vector<int>::iterator it =
	       arrangementFallbackOptions[j].begin();
	     it != arrangementFallbackOptions[j].end(); ++it){
	  if (*it == i){
	    expr += myVars[j + 1];
	  }
	}
      }
      std::string s = "vertex_" + KPDFunctions::intToString(i);
      model.addConstr(expr, GRB_LESS_EQUAL, 1.0, s);
    }

    model.optimize();

    int optimstatus = model.get(GRB_IntAttr_Status);
    double objval = 0;
    if (optimstatus == GRB_OPTIMAL) { // GRB_OPTIMAL when the MIPGap is reached, default is 1e-4;
      objval = model.get(GRB_DoubleAttr_ObjVal);
    }
    else if (optimstatus == GRB_INF_OR_UNBD) {
      std::cout << "Model is infeasible or unbounded" << std::endl; return;
    }
    else if (optimstatus == GRB_INFEASIBLE) {
      std::cout << "Model is infeasible" << std::endl; return;
    }
    else if (optimstatus == GRB_UNBOUNDED) {
      std::cout << "Model is unbounded" << std::endl; return;
    }
    else {
      std::cout << "Optimization was stopped with status = " << optimstatus << std::endl; return;
    }

    int numvars = model.get(GRB_IntAttr_NumVars);
    vars = model.getVars();

    for (int j = 0; j < numvars; j++) {
      GRBVar v = vars[j];
      if (v.get(GRB_DoubleAttr_X) == 1){ selectedFallbacks.push_back(j); }
    }

  }
  catch (GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
  }
  catch (...) {
    std::cout << "Error during optimization" << std::endl;
  }

  delete[] myVars;
  delete[] vars;
  delete env;

  for (std::vector<int>::iterator fallbackOptionIt = selectedFallbacks.begin();
       fallbackOptionIt != selectedFallbacks.end(); fallbackOptionIt++) {
    std::vector<int> cycleOrChain;
    for (std::vector<int>::iterator nodeIt =
	   arrangementFallbackOptions[*fallbackOptionIt].begin();
	 nodeIt != arrangementFallbackOptions[*fallbackOptionIt].end(); nodeIt++) {
      cycleOrChain.push_back(arrangementNodes[*nodeIt]->getID());
    }
    options.push_back(cycleOrChain);
  }

  
  // Sanity check fallback solution
  if (!KPDFunctions::isValidSolution(options))
    throw std::logic_error ("Invalid fallbacks. Simulation terminated\n");

}  // 












std::string KPDArrangement::toLabel(const std::string &_delim) {
  std::stringstream ss;
  for (std::vector<KPDNode *>::iterator it = arrangementNodes.begin();
       it != arrangementNodes.end() - 1; it++)
    ss << (*it)->getID() << _delim;
  
  ss << arrangementNodes.back()->getID();
  return ss.str();
}






// KPDArrangement::toTransplantString()
// -------------------------------------------------------------------
std::string KPDArrangement::toTransplantString() {
  std::stringstream ss;
  int N = (int)arrangementNodes.size();
  int adjacencyCount = 0;
  
  for (int i = 1; i <= N; i++) {
    int donorNodeIndex = i - 1;

    KPDNode * donorNode = arrangementNodes[donorNodeIndex];
    KPDNodeType donorNodeType = arrangementNodeTypes[donorNodeIndex];

    int numDonors = donorNode->getNumberOfDonors();

    for (int j = 1; j <= N; j++) {
      int candidateNodeIndex = j - 1;

      if (arrangementNodeTypes[candidateNodeIndex] == PAIR &&
	  arrangementMatches.find(donorNodeIndex) != arrangementMatches.end() &&
	  arrangementMatches[donorNodeIndex].find(candidateNodeIndex) !=
	  arrangementMatches[donorNodeIndex].end() &&
	  !arrangementMatches[donorNodeIndex][candidateNodeIndex].empty()
	  ) {
	KPDNode * candidateNode = arrangementNodes[candidateNodeIndex];

	for (int donorIndex = 0; donorIndex < numDonors; donorIndex++) {
	  KPDMatch * match =
	    arrangementMatches[donorNodeIndex][candidateNodeIndex][donorIndex]->copy();

	  adjacencyCount += match->getAdjacency();

	  if (match->getAdjacency()) {
	    ss << currentIteration << ","
	       << matchRunTime << ","
	       << transplantationTime << ","
	       << KPDFunctions::nodeTypeToString(donorNodeType) << ","
	       << donorNode->getID() << ","
	       << (donorIndex + 1) << ","
	       << donorNode->getArrivalTime() << ","
	       << donorNode->getDonorString(donorIndex) << ","
	       << candidateNode->getID() << ","
	       << candidateNode->getArrivalTime() << ","
	       << candidateNode->getNumberOfDonors() << ","
	       << KPDFunctions::boolToYesNo(arrangementAvailability[j]) << ","
	       << candidateNode->getCandidateString() << ","
	       << match->matchString() << ","
	       << KPDFunctions::boolToYesNo(arrangementTransplantSelected
					    [donorNodeIndex][candidateNodeIndex] == donorIndex)
	       << ","
	       << _tempApproxPoolSize
	       << "\n";
	  }
	}
      }  // if (arrangementNodeTypes[candidateNodeIndex] == PAIR)
    }  // for (int j = 1; j <= N; j++)
  }  // for (int i = 1; i <= N; i++)

  // std::cout << "# of adjacencies (DD-Arrangement): " << adjacencyCount << std::endl;

  return ss.str();
}


#endif
