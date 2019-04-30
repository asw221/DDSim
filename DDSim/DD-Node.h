/* ---------------------------------------------
DD-Node.h
Forms KPD nodes (donor-candidate pairings or NDDs)
---------------------------------------------- */

#ifndef NODE_H
#define NODE_H

#include "DD-Enums-Functions.h"
#include "DD-Candidate.h"
#include "DD-Donor.h"

#include <string>
#include <vector>

class KPDNode {
	
private:

  int nodeID;

  std::vector<KPDDonor *> nodeDonors;
  KPDCandidate * nodeCandidate;

  int nodeArrivalTime;
  KPDNodeType nodeType;
	
public:
	
  //Constructors

  //Construct a dummy node
  KPDNode();
  //Construct an NDD
  KPDNode(int id, int arrivalTime, KPDDonor * donor);	
  //Construct a donor-candidate pair
  KPDNode(int id, int arrivalTime, std::vector<KPDDonor * > donors, KPDCandidate * candidate);
  ~KPDNode();

  //IDs
  int getID();
  void setID(int id);
	
  //Candidate Information
  KPDCandidate * getCandidate();

  KPDBloodType getCandidateBT();
  int	getCandidatePRA();

  void setCandidate(KPDCandidate * candidate);

  //Donor Information
  KPDDonor * getDonor(int donorIndex);
  std::vector<KPDDonor *> getDonors();
  int getNumberOfDonors();

  KPDBloodType getDonorBT(int donorIndex);	

  void setDonors(std::vector<KPDDonor *> donors);
	
  //Node Type Information
  KPDNodeType getType();
  void setType(KPDNodeType type);
  bool isNDD() const;
  bool isDeceasedDonor() const;

  //Node Arrival Information
  int getArrivalTime();		
  void setArrivalTime(int time);

  //Strings
  std::string getCandidateString(const bool &header = false);
  std::string getDonorString(const int &donorIndex, const bool &header = false);

  //Copy Constructor
  KPDNode * copy();
};




KPDNode::KPDNode() { ; }




KPDNode::KPDNode(int id, int arrivalTime, KPDDonor * donor) {
  nodeID = id;
  nodeDonors.push_back(donor->copy());
  nodeArrivalTime = arrivalTime;
  nodeType = NDD;
}




KPDNode::KPDNode(
  int id,
  int arrivalTime,
  std::vector<KPDDonor *> donors,
  KPDCandidate * candidate
) {
  nodeID = id;
  for (std::vector<KPDDonor *>::iterator it = donors.begin();
       it != donors.end(); it++) {
    nodeDonors.push_back((*it)->copy());
  }
  nodeCandidate = candidate->copy();
  nodeArrivalTime = arrivalTime;
  nodeType = PAIR;
}




KPDNode::~KPDNode(){
  nodeDonors.clear();
  delete nodeCandidate;
}




int KPDNode::getID() {
  return nodeID;
}

void KPDNode::setID(int id){
  nodeID = id;
}

KPDCandidate * KPDNode::getCandidate() {
  return nodeCandidate->copy();
}

KPDBloodType KPDNode::getCandidateBT() {
  return nodeCandidate->getBT();
}

int	KPDNode::getCandidatePRA() {
  return nodeCandidate->getPRA();
}

void KPDNode::setCandidate(KPDCandidate * candidate) {
  nodeCandidate = candidate;
}

KPDDonor * KPDNode::getDonor(int donorIndex) {	
  return nodeDonors.at(donorIndex)->copy();
}




std::vector<KPDDonor *> KPDNode::getDonors() {
  std::vector<KPDDonor *> donors;
  for (std::vector<KPDDonor *>::iterator it = nodeDonors.begin(); it != nodeDonors.end(); it++) {
    donors.push_back((*it)->copy());
  }

  return donors;
}




int KPDNode::getNumberOfDonors() {
  return (int)nodeDonors.size();
}




KPDBloodType KPDNode::getDonorBT(int donorIndex) {
  return nodeDonors.at(donorIndex)->getBT();
}




void KPDNode::setDonors(std::vector<KPDDonor *> donors) {
  nodeDonors.clear();
  for (std::vector<KPDDonor *>::iterator it = donors.begin(); it != donors.end(); it++) {
    nodeDonors.push_back(*it);
  }
}




KPDNodeType KPDNode::getType() {
  return nodeType;
}




void KPDNode::setType(KPDNodeType type) {
  nodeType = type;
}




bool KPDNode::isNDD() const {
  return nodeType != PAIR;
};


bool KPDNode::isDeceasedDonor() const {
  bool isDeceased = false;
  // std::vector<KPDDonor *> nodeDonors;
  for (int i = 0; i < nodeDonors.size(); i++)
    isDeceased = isDeceased || nodeDonors[i]->isDeceasedDonor();
  return isDeceased;
};


int KPDNode::getArrivalTime(){
  return nodeArrivalTime;
}

void KPDNode::setArrivalTime(int time){
  nodeArrivalTime = time;
}




std::string KPDNode::getCandidateString(const bool &header) {
  if (nodeType == PAIR) {
    return nodeCandidate->candidateString(header);
  }
  else {
    return ",,,,,,,,,,,,"; // No Candidate Info for NDD
    // return KPDFunctions::csvStringToBlank(nodeCandidate->candidateString(true));
  }
}




std::string KPDNode::getDonorString(const int &donorIndex, const bool &header) {	
  return nodeDonors[donorIndex]->donorString(header);
}




KPDNode * KPDNode::copy() {
  KPDNode * copyNode;
  if (nodeType == NDD){
    copyNode = new KPDNode(nodeID, nodeArrivalTime, nodeDonors.front()->copy());
  }
  else {
    std::vector<KPDDonor *> copiedDonors;
    for (std::vector<KPDDonor*>::iterator it = nodeDonors.begin(); it != nodeDonors.end(); it++) {
      copiedDonors.push_back((*it)->copy());
    }
    copyNode = new KPDNode(nodeID, nodeArrivalTime, copiedDonors, nodeCandidate->copy());
  }

  return copyNode;
}

#endif
