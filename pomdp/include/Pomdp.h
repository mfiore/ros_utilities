/*
 * Pomdp.h
 *
 *  Created on: May 9, 2015
 *      Author: mfiore
 */

#ifndef SOURCE_DIRECTORY__UTILITIES_POMDP_POMDP_H_
#define SOURCE_DIRECTORY__UTILITIES_POMDP_POMDP_H_

/*
 * File:   Pomdp.h
 * Author: mfiore
 *
 * Created on August 18, 2014, 11:53 AM
 */

#ifndef POMDP_H
#define	POMDP_H

#include <vector>
#include<map>

#include "VariableSet.h"
#include <fstream>
#include "PriorityQueue.h"
#include <tuple>

#include "NestedLoop.h"

#include "AlphaVector.h"
using namespace std;



class Pomdp {
public:
    Pomdp();
    Pomdp(const Pomdp& orig);

    virtual ~Pomdp();

    void createPomdp(string name, bool rewrite=false);
    void setInitialState(std::vector<pair<VariableSet, double>> initialBelief);

    //three functions to override in the derived classes
    virtual std::map<VariableSet, double> transitionFunction(VariableSet state, string action) = 0;
    virtual double observationFunction(VariableSet observationState, VariableSet state, string action)=0;

    virtual int rewardFunction(VariableSet state, string action) = 0;

    map<int,double> updateBelief(map<int,double> belief,string action, map<string,string> obsStates,
    		int o);
    void updateCurrentBelief(string action, std::map<string, string> obsStates, std::map<string, string> observations);

    map<int,double> filterBelief(map<int,double> belief,map<string,string> obsStates);

    //learning functions
    double updateValue(std::map<int,double> belief);
    void updateAlphaVectors();
    void valueIteration(bool rewrite=false);
    vector<map<int,double>> expandBeliefSet(vector<map<int,double>> beliefSet);
	vector<AlphaVector> improveValue(vector<AlphaVector> vAlpha,vector<map<int,double>> bSet);
	AlphaVector backup(map<int,double> b,vector<AlphaVector> vAlpha,
			std::map<std::tuple<int,string,int,int>,double> alpha_aos);
	std::map<std::tuple<int,string,int,int>,double> calculateAlphaAOS(int i, AlphaVector alpha,
			std::map<std::tuple<int,string,int,int>,double> currentSet);


    string chooseAction();

    //utility functions
    void printBelief();
    void printTransitionFunction();
    void printRewardFunction();
    void printStates();


    //podmp specification
    //includes all the variables of the system (observed state variables, hidden state variables, observations
    //and the action var) with their values
    std::vector<string> observedVariables_;
    std::vector<string> hiddenVariables_;
    std::vector<string> observationsVariables_;

    std::vector<string> variables_;

    std::map<string, std::vector<string>> varValues_;
    std::vector<string> actions_;

    std::map<string,std::vector<string>> observationValues_;


    //variables enumerations
    //the system state enumeration is kept both way (from enumeration number to system state and other way)
    std::map<VariableSet, int> mapStateEnum_;
    std::vector<VariableSet> vecStateEnum_;

    std::map<VariableSet,int> mapObservationEnum_;
    std::vector<VariableSet>  vecObsEnum_;

    std::map<pair<int, string>, std::map<int, double>> transitionMap_;

    std::map<pair<int, string>, int> rewardMap_; //reward function
    std::map<std::tuple<int,int,string>, double> observationMap_;

    std::map<int, double> currentBelief_;


    std::vector<AlphaVector> alphaVectors_;

    string name_;

private:


};

#endif	/* POMDP_H */
#endif /* SOURCE_DIRECTORY__UTILITIES_POMDP_POMDP_H_ */
