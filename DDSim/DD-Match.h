/* ---------------------------------------------
DD-Match.h
Stores information on donor-candidate matches
(within KPD, DD to KPD, and KPD to Waitlist)
---------------------------------------------- */

#ifndef MATCH_H
#define MATCH_H

#include "DD-Enums-Functions.h"

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

class KPDMatch {

private:
  bool adjacency; // This is the assumed compatibility of the match given no other information

  double fiveYearSurvival;
  double tenYearSurvival;
  double transplantDifficultyScore;
  double randomUtility;

  double assumedSuccessProbability;
  double actualSuccessProbability;

  KPDCrossmatch virtualCrossmatchResult;
	
  bool successfulMatch; // This is the true underlying compatibility of the match

public:
  //Constructors
  //Construts a dummy match
  KPDMatch();
  
  //Constructs a match with all defined match characteristics 
  KPDMatch(
    bool match,
    double fiveYearSurv,
    double tenYearSurv,
    double score,
    double util,
    double assumedProb,
    double actualProb,
    KPDCrossmatch virtualResult,
    bool success
  );
  
  ~KPDMatch();

  //Getters
  bool getAdjacency();

  double getFiveYearSurvival();
  double getTenYearSurvival();
  double getTransplantDifficultyScore();
  double getRandomUtility();

  double getUtility(KPDUtilityScheme scheme);
  // Returns utility value based on specified utility scheme

  double getAssumedSuccessProbability();
  double getActualSuccessProbability();

  KPDCrossmatch getVirtualCrossmatchResult();
  bool getSuccessfulMatch();


  //Setters
  void setAdjacency(bool match);

  void setFiveYearSurvival(double surv);
  void setTenYearSurvival(double surv);
  void setTransplantDifficultyScore(double score);
  void setRandomUtility(double util);

  void setUtility(KPDUtilityScheme scheme, double util);
  // Sets utility value based on specified utility scheme

  void setAssumedSuccessProbability(double prob);
  void setActualSuccessProbability(double prob);

  void setVirtualCrossmatchResult(KPDCrossmatch result);
  void setSuccessfulMatch(bool success);
	
  void setMatchProperties(
    bool match,
    double fiveYearSurv,
    double tenYearSurv,
    double score,
    double util,
    double assumedProb,
    double actualProb,
    KPDCrossmatch virtualResult,
    bool success
  );

  //Strings
  std::string matchString(); // Returns a comma-separated string
  std::string matchOutput(); // Returns detailed match information
  std::string matchShortOutput(); // Returns brief match information

  //Copy Constructor
  KPDMatch * copy() const;
  
};




KPDMatch::KPDMatch() {
  adjacency = false;

  fiveYearSurvival = 0.0;
  tenYearSurvival = 0.0; 
  transplantDifficultyScore = 0.0;
  randomUtility = 0.0;
	
  assumedSuccessProbability = 0.0;
  actualSuccessProbability = 0.0;

  virtualCrossmatchResult = CROSSMATCH_FAILED;

  successfulMatch = false;
}




KPDMatch::KPDMatch(
  bool match,
  double fiveYearSurv,
  double tenYearSurv,
  double score,
  double util,
  double assumedProb,
  double actualProb,
  KPDCrossmatch virtualResult,
  bool success
) {	
  adjacency = match;

  fiveYearSurvival = fiveYearSurv;
  tenYearSurvival = tenYearSurv;
  transplantDifficultyScore = score;
  randomUtility = util;

  assumedSuccessProbability = assumedProb;
  actualSuccessProbability = actualProb;

  virtualCrossmatchResult = virtualResult;

  successfulMatch = success;
};




KPDMatch::~KPDMatch() { ; }

bool KPDMatch::getAdjacency(){
  return adjacency;
}

double KPDMatch::getFiveYearSurvival(){
  return fiveYearSurvival;
}

double KPDMatch::getTenYearSurvival(){
  return tenYearSurvival;
}

double KPDMatch::getTransplantDifficultyScore(){
  return transplantDifficultyScore;
}

double KPDMatch::getRandomUtility(){
  return randomUtility;
}

double KPDMatch::getUtility(KPDUtilityScheme scheme){
  if (scheme == UTILITY_FIVE_YEAR_SURVIVAL) {
    return fiveYearSurvival;
  }
  else if (scheme == UTILITY_TEN_YEAR_SURVIVAL) {
    return tenYearSurvival;
  }
  else if (scheme == UTILITY_TRANSPLANT_DIFFICULTY) {
    return transplantDifficultyScore;
  }
  else if (scheme == UTILITY_RANDOM) {
    return randomUtility;
  }
  else {
    return 1.0;
  }
}

double KPDMatch::getAssumedSuccessProbability(){
  return assumedSuccessProbability;
}

double KPDMatch::getActualSuccessProbability(){
  return actualSuccessProbability;
}

KPDCrossmatch KPDMatch::getVirtualCrossmatchResult(){
  return virtualCrossmatchResult;
}

bool KPDMatch::getSuccessfulMatch(){
  return successfulMatch;
}




void KPDMatch::setAdjacency(bool match){
  adjacency = match;
}

void KPDMatch::setFiveYearSurvival(double surv){
  fiveYearSurvival = surv;
}

void KPDMatch::setTenYearSurvival(double surv){
  tenYearSurvival = surv;
}

void KPDMatch::setTransplantDifficultyScore(double score){
  transplantDifficultyScore = score;
}

void KPDMatch::setRandomUtility(double util){
  randomUtility = util;
}

void KPDMatch::setUtility(KPDUtilityScheme scheme, double util){
  if (scheme == UTILITY_FIVE_YEAR_SURVIVAL) {
    fiveYearSurvival = util;
  }
  else if (scheme == UTILITY_TEN_YEAR_SURVIVAL) {
    tenYearSurvival = util;
  }
  else if (scheme == UTILITY_TRANSPLANT_DIFFICULTY) {
    transplantDifficultyScore = util;
  }
  else if (scheme == UTILITY_RANDOM) {
    randomUtility = util;
  }
}

void KPDMatch::setAssumedSuccessProbability(double prob){
  assumedSuccessProbability = prob;
}

void KPDMatch::setActualSuccessProbability(double prob){
  actualSuccessProbability = prob;
}

void KPDMatch::setVirtualCrossmatchResult(KPDCrossmatch result){
  virtualCrossmatchResult = result;
}

void KPDMatch::setSuccessfulMatch(bool success){
  successfulMatch = success;
}




void KPDMatch::setMatchProperties(
  bool match, double fiveYearSurv, double tenYearSurv,
  double score, double util,
  double assumedProb, double actualProb,
  KPDCrossmatch virtualResult, bool success
) {	
  adjacency = match;
  fiveYearSurvival = fiveYearSurv;
  tenYearSurvival = tenYearSurv;
  transplantDifficultyScore = score;
  randomUtility = util;
  assumedSuccessProbability = assumedProb;
  actualSuccessProbability = actualProb;
  virtualCrossmatchResult = virtualResult;
  successfulMatch = success;
}




std::string KPDMatch::matchString(){
  std::stringstream matchInfo;
  matchInfo << fiveYearSurvival << ",";
  matchInfo << tenYearSurvival << ",";						
  matchInfo << transplantDifficultyScore << ",";
  matchInfo << randomUtility << ",";
  matchInfo << assumedSuccessProbability << ",";
  matchInfo << actualSuccessProbability << ",";
  matchInfo << KPDFunctions::crossmatchToString(virtualCrossmatchResult) << ",";
  matchInfo << KPDFunctions::boolToTF(successfulMatch);

  return matchInfo.str();
}




std::string KPDMatch::matchOutput(){
  int tab = 3;
  std::stringstream matchInfo;
  matchInfo << KPDFunctions::indent(tab)
	    << "5-Year Survival: " << fiveYearSurvival
	    << ", 10-Year Survival: " << tenYearSurvival
	    << ", Score: " << transplantDifficultyScore
	    << ", Random Utility: "
	    << randomUtility
	    << std::endl;
  
  matchInfo << KPDFunctions::indent(tab)
	    << "Assumed Success Probability: " << assumedSuccessProbability
	    << ", Actual Success Probability: " << actualSuccessProbability
	    << std::endl;
  
  matchInfo << KPDFunctions::indent(tab)
	    << "Virtual Crossmatch Result: "
	    << KPDFunctions::crossmatchToString(virtualCrossmatchResult)
	    << std::endl;
  
  matchInfo << KPDFunctions::indent(tab)
	    << "Successful Match (Hidden): " << KPDFunctions::boolToTF(successfulMatch)
	    << std::endl;

  return matchInfo.str();
}




std::string KPDMatch::matchShortOutput() {
  std::stringstream matchInfo;
  matchInfo << "["
	    << fiveYearSurvival << ","
	    << tenYearSurvival << ","
	    << transplantDifficultyScore << ","
	    << randomUtility << ",";
  matchInfo << assumedSuccessProbability << ","
	    << actualSuccessProbability << ",";
  matchInfo << KPDFunctions::crossmatchToString(virtualCrossmatchResult) << ","
	    << KPDFunctions::boolToTF(successfulMatch) << "]";

  return matchInfo.str();
}




KPDMatch * KPDMatch::copy() const {
  KPDMatch * copyMatch = new KPDMatch(
    adjacency, fiveYearSurvival, tenYearSurvival,
    transplantDifficultyScore, randomUtility,
    assumedSuccessProbability, actualSuccessProbability,
    virtualCrossmatchResult, successfulMatch
  );

  return copyMatch;
}


#endif
