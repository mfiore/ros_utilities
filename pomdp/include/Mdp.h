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



class Mdp {
public:
    Mdp();
    Mdp(const Mdp& orig);
    //    Pomdp(std::vector <PomdpVar> obsStates, std::vector<PomdpVar> hiddenStates, std::vector<string> actions, std::vector<PomdpVar> observations,
    //    std::map<string, std::vector< pair<PomdpVar, double > > > transition, std::map <pair<string,string>, double> observationProb, std::map<string, int> reward, 
    //    double discount, std::map<std::vector<PomdpVar> , double > initialBelief);
    virtual ~Mdp();

    void createMdp(string name, bool rewrite=false);
    void setInitialState(std::vector<pair<VariableSet, double>> initialBelief);

    //three functions to override in the derived classes
    virtual std::map<VariableSet, double> transitionFunction(VariableSet state, string action) = 0;
    virtual int rewardFunction(VariableSet state, string action) = 0;
    virtual bool isGoalState(VariableSet state) = 0;


    map<int,double> pomdpUpdateBelief(map<int,double> belief,string action, std::map<string, string> obsStates,
    		std::map<string, string> observations);
    void pomdpUpdateCurrentBelief(string action, std::map<string, string> obsStates, std::map<string, string> observations);
    void updateCurrentBelief(string action, std::map<string, string> obsStates);

    map<int,double> filterBelief(map<int,double> belief,std::map<string, string> observedVariables);

    std::map<string,double> getNormQValue();

    double getQValue(string action);
    
    //learning functions
    int bellmanBackup(int i, std::vector<double> vhi);
    void valueIteration(bool rewrite=false);
    void prioritizedSweeping(); //for now it's slower on the tested examples. Maybe it's the queue overhead
    double pomdpUpdateValue(std::map<int,double> belief);
    void pomdpUpdateAlphaVectors();

    void pomdpValueIteration(bool rewrite=false);
    vector<map<int,double>> expandBeliefSet(vector<map<int,double>> beliefSet);
	vector<AlphaVector> improveValue(vector<AlphaVector> vAlpha,vector<map<int,double>> bSet);
	AlphaVector backup(map<int,double> b,vector<AlphaVector> vAlpha,
			std::map<std::tuple<int,string,int,int>,double> alpha_aos);
	std::map<std::tuple<int,string,int,int>,double> calculateAlphaAOS(int i, AlphaVector alpha,
			std::map<std::tuple<int,string,int,int>,double> currentSet);


    string chooseAction();
    double getActionValue(string action);

    //utility functions
    void printBelief();
    void printTransitionFunction();
    void printRewardFunction();
    void printStates();
    void printActualQValues();



    //Variables
    std::map<VariableSet, int> mapStateEnum; //the system state enumeration is kept both way (from enumeration number to system state and other way)
    std::vector<VariableSet> vecStateEnum;
    std::map<VariableSet,int> mapObsEnum;
    std::vector<VariableSet>  vecObsEnum;

    std::map<pair<int, string>, std::map<int, double>> transition;
    std::map<pair<int, string>, std::vector<int>> predecessors; //enumeration of the transition function
    std::map<pair<int, string>, int> reward; //reward function
    std::map<std::tuple<int,int,string>, double> observationFunction;
    std::vector<int> goalStates;


    //podmp specification
    //includes all the variables of the system (observed state variables, hidden state variables, observations
    //and the action var) with their values
    std::map<string, std::vector<string> > varValues;
    std::vector<string> variables; //lists all variables;
    std::vector<string> observedVariables; //lists which variables are observed

    std::map<string, std::vector<string>> observationValues;
    std::vector<string> observationVariables; //lists which variables are observations
    std::vector<string> hiddenVariables; //lists which variables are hidden
    std::vector<string> actions;


    std::map<pair<int, string>, double> qValue; //human action values
    std::map<int, double> belief;


    std::vector<std::vector<double> > alphaVectors;

    string name;

private:


};

#endif	/* POMDP_H */

