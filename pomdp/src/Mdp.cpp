/* 
 * File:   Mdp.cpp
 * Author: mfiore
 * 
 * Created on August 18, 2014, 11:53 AM
 */

#include "../include/Mdp.h"

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <stack>
#include <utility>
#include <algorithm>
#include <c++/4.6/bits/stl_queue.h>

Mdp::Mdp() {

    //    this->varValues=varValues;
    //    for (auto var:varValues) {
    //        variables.push_back(var.first);
    //    }
    //    this->hiddenVariables=hiddenVariables;
    //    this->observedVariables=observedVariables;
    //    for (auto var:observationValues) {
    //        observationVariables.push_back(var.first);
    //    }
    //    this->observationValues=observationValues;
    //    this->actions=actions;
}

Mdp::Mdp(const Mdp& orig) {
}

Mdp::~Mdp() {
}



/* This function filters the current belief with the observed states, erasing the belief that don't respect them.
 
 * obsStates: the observed states and their values
 */
map<int,double> Mdp::filterBelief(map<int,double> belief,std::map<string, string> obsStates) {
    std::map<int, double> newBelief;
    bool checkObservedVariables = true;

    for (auto currentBelief : belief) {
        checkObservedVariables = true;

        //we check if all the observed variables in the belief are equal the observed variables in input
        VariableSet v = vecStateEnum[currentBelief.first];
        for (auto j : obsStates) {
            if (v.set.at(j.first) != j.second) {
                checkObservedVariables = false;
                break;
            }
        }
        //if so we keep the belief
        if (checkObservedVariables) {
            newBelief[currentBelief.first] = currentBelief.second;

        }
    }
    //normalize
    double sum = 0;
    for (auto j : newBelief) {
        sum = sum + j.second;
    }
    for (auto j : newBelief) {
        j.second = j.second / sum;
    }

    return newBelief;

}


map<int,double> Mdp::pomdpUpdateBelief(map<int,double> belief,string action, std::map<string, string> obsStates, std::map<string, string> observations) {
	  std::map<int, double> finalBelief;

	    VariableSet vso;
	    vso.set=observations;
	    int o=mapObsEnum[vso];

	    //for each current belief
	    for (auto beliefState: belief) {
	        int s = beliefState.first;

	        pair<int, string> transitionInput{s, action};
	        std::map<int, double> futureBelief = transition[transitionInput];

	        for (auto futureBeliefState: futureBelief) {
	            finalBelief[futureBeliefState.first] = finalBelief[futureBeliefState.first] + futureBeliefState.second;
	        }

	        //normalize
	        double sum = 0;
	        for (auto aBelief : finalBelief) {
	            sum = sum + aBelief.second;
	        }
	        for (auto aBelief : finalBelief) {
	            aBelief.second = aBelief.second / sum;
	        }
	    }
	    //apply observation function
	    for (auto beliefState:finalBelief) {
	    	std::tuple<int,int,string> observationInput{o,beliefState.first,action};
	    	beliefState.second=beliefState.second*observationFunction[observationInput];
	    }
	    //normalize
	     double sum = 0;
	     for (auto aBelief : finalBelief) {
	         sum = sum + aBelief.second;
	     }
	     for (auto aBelief : finalBelief) {
	         aBelief.second = aBelief.second / sum;
	     }
	  	finalBelief=filterBelief(finalBelief,obsStates);

	    return finalBelief;
}


/*
 This function updates the current system belief for a POMDP
 * 
 * action: current action
 * obsStates: set of observed variables
 * observations: set of observations
 */
void Mdp::pomdpUpdateCurrentBelief(string action, std::map<string, string> obsStates, std::map<string, string> observations) {
	map<int,double> newBelief;
	newBelief=pomdpUpdateBelief(this->belief,action,obsStates,observations);
    this->belief=newBelief;
}


/*
 This function updates the current system belief for an MDP
 *
 * action: current action
 * obsStates: set of observed variables
 */
void Mdp::updateCurrentBelief(string action, std::map<string, string> obsStates) {
    std::map<int, double> finalBelief;

    //for each current belief
    for (auto currentBelief : belief) {
        int s = currentBelief.first;

        pair<int, string> transitionInput{s, action};
        std::map<int, double> futureBeliefs = transition[transitionInput];

        for (auto aBelief : futureBeliefs) {
            finalBelief[aBelief.first] = finalBelief[aBelief.first] + aBelief.second;
        }
        //normalize
        double sum = 0;
        for (auto aBelief : finalBelief) {
            sum = sum + aBelief.second;
        }
        for (auto aBelief : finalBelief) {
            aBelief.second = aBelief.second / sum;
        }
    }
    finalBelief=filterBelief(finalBelief,obsStates);
    this->belief = finalBelief;


}


/*
 Executes a bellman backup on state i that is max a [R(i,a) + L * sum T(i,a,s) vhi(s)]
 * returns the new v values and updates the q values
 */
int Mdp::bellmanBackup(int i, std::vector<double> vhi) {
    double maxValue = 0;
    for (string action : actions) { //for every action
        pair<int, string> rInput{i, action}; //calculate the reward of this state with this action

        int currentReward;
        //we have to do this because the std::map increments of one element if we just use [] with a non existing member and the reward function
        //contains only transition with non zero rewards.
        if (reward.find(rInput) != reward.end()) {
            currentReward = reward[rInput];
        } else {
            currentReward = 0;
        }
        pair<int, string> transitionInput{i, action};
        std::map<int, double> futureBeliefs = transition[transitionInput];

        double sum = 0;
        for (auto aBelief : futureBeliefs) {
            sum = sum + aBelief.second * vhi[aBelief.first]; //sum on the probabilities of the future states * the value of reaching that state
        }
        pair<int, string> qInput{i, action};
        double havNew = currentReward + 0.3 * sum; //0.3 weights the future rewards
        qValue[qInput] = havNew; //update the human action value
        if (qValue[qInput] > maxValue) {
            maxValue = qValue[qInput];
        }

    }
    return maxValue;
}

/*
 This procedure does value iteration and computes the Q-Values for the actions. At the moment we don't use the policy learnt through here,
 * preferring to use the one with sarsop. The algorithm works on enumerations of states
 */
void Mdp::valueIteration(bool rewrite) {
    string fileName = name + ".policy";
    ifstream inputFile(fileName);
    if (inputFile.good() && !rewrite) {
        string line;
        for (int i = 0; i < vecStateEnum.size(); i++) {
            int i1 = 0;
            getline(inputFile, line);
            for (string action : actions) {
                int i2 = line.find(" ", i1);
                string qv = line.substr(i1, i2 - i1);
                pair<int, string> qInput{i, action};
                qValue[qInput] = stod(qv);
                i1 = i2 + 1;
            }
        }
        inputFile.close();

    } else {
        double epsilon = 0.1; //stopping value of the value iteration
        std::vector<double> vhi(vecStateEnum.size(), 0); //vhi mantains the predicted reward in a state following the optimal policy
        std::vector<double> vhiOld(vecStateEnum.size(), 0);
        double maxDiff = 0;
        cout << "Starting Value Iteration\n";
        do { //we loop until vhiOld-hvi<epsilon (more or less)

            for (int s = 0; s < vecStateEnum.size(); s++) { //for each enumeration
                vhiOld[s] = vhi[s];
                if (s == 2) {
                    cout << "ah";
                }
                vhi[s] = bellmanBackup(s, vhi);
                //            if (vhi[s]-vhiOld[s]==7) {
                //                cout<<"ah";
                //            }
            }
            maxDiff = 0;
            for (int i = 0; i < vhi.size(); i++) { //calculate the maximum difference on the vhi (stopping parameter)
                double d = abs(vhi[i] - vhiOld[i]);
                if (d > maxDiff) {
                    maxDiff = d;
                    cout << "Difference " << d << "\n";
                }
            }

        } while (maxDiff > epsilon);
        ofstream file(fileName);
        for (int i = 0; i < vecStateEnum.size(); i++) {
            for (string action : actions) {
                pair<int, string> qInput{i, action};
                file << qValue[qInput] << " ";
            }
            file << "\n";
        }
        file.close();
    }
}

void Mdp::prioritizedSweeping() {
    std::vector<double> vhi(vecStateEnum.size(), 0);
    std::vector<double> vhiOld(vecStateEnum.size(), 0);
    double epsilon = 0.1;
    PriorityQueue pq;
    for (int i = 0; i < vecStateEnum.size(); i++) {
        pq.pushElement(pair<int, double>{i, 1});
    }
    //    for (int i = 0; i < goalStates.size(); i++) {
    //        pq.pushElement(pair<int, double>{i, 1});
    //    }

    while (!pq.isEmpty()) {
        pair<int, double> actualElement = pq.pop();
        int s = actualElement.first;
        //        cout << pq.size() << "\n";
        //
        //        cout << s << "\n";
        vhiOld[s] = vhi[s];
        vhi[s] = bellmanBackup(s, vhi);

        std::map<int, bool> isVisited;
        for (string action : actions) {
            pair<int, string> predInput{s, action};
            std::vector<int> vsp = predecessors[predInput];

            for (int sp : vsp) {
                if (isVisited[sp] == false) {
                    isVisited[sp] = true;

                    double maxValue = 0;
                    for (string action2 : actions) {
                        pair<int, string> transitionInput{sp, action2};
                        std::map<int, double> futureStates = transition[transitionInput];
                        double value = futureStates[s] * (vhi[s] - vhiOld[s]);
                        if (value > maxValue) {
                            maxValue = value;
                        }
                    }
                    if (sp != s) {
                        double sppriority = pq.getPriority(sp);
                        if (sppriority > maxValue) {
                            maxValue = sppriority;
                        }
                    }
                    if (maxValue > epsilon) {
                        //                        cout << maxValue << "\n";

                        pair<int, double> newElement{sp, maxValue};
                        pq.pushElement(newElement);
                        //                        pq.print();
                    }
                }

            }
        }
    }
}


std::map<std::tuple<int,string,int,int>,double> Mdp::calculateAlphaAOS(int i, AlphaVector alpha,
		std::map<std::tuple<int,string,int,int>,double> currentSet) {
	for (string a:actions) {
		for (int o=0; o<vecObsEnum.size(); o++) {
			for (int s=0; s<vecStateEnum.size(); s++) {
				double sum=0;
				for (int sp=0;s<vecStateEnum.size();s++) {
					std::tuple<int,int,string> observationInput=std::make_tuple(o,sp,a);
					pair<int,string> transitionInput{s,a};
					sum=sum+alpha.coefficients[sp]*observationFunction[observationInput]*transition[transitionInput][sp];
				}
				auto alphaIndex=std::make_tuple(i,a,o,s);
				currentSet[alphaIndex]=sum;
			}
		}
	}
	return currentSet;
}

AlphaVector Mdp::backup(map<int,double> b,vector<AlphaVector> vAlpha,
		std::map<std::tuple<int,string,int,int>,double> alpha_aos) {

	double gamma=0.3;
	//compute alpha_ab
	vector<AlphaVector> v_alpha_ab;
	//there is a different alpha_ab for each alpha and for each a
	for (string a:actions) {
		  for (AlphaVector alpha:vAlpha) {
			  AlphaVector sumO;
			  sumO.action=a;
			  for (int o=0; o<vecObsEnum.size();o++) {
				  //find argmax a_o| a_o*b
				  double max_alpha_ao;
				  int max_alpha_index=0;
				  for (int alphaIndex=0; alphaIndex<vAlpha.size();alphaIndex++) {
					  //get this alpha_ao
					  double sum=0;
					  for (int s=0; s<vecStateEnum.size();s++) {
						  auto alpha_ao_index=std::make_tuple(alphaIndex,a,o,s);
						  sum=b[s]*alpha_aos[alpha_ao_index];
					  }
					  if (sum>max_alpha_ao) {
						  max_alpha_ao=sum;
						  max_alpha_index=alphaIndex;
					  }
				  }
				  //sum argMax to current max sum
				  for (int s=0; s<vecStateEnum.size();s++) {
					  auto alpha_ao_index=std::make_tuple(max_alpha_index,a,o,s);
					  sumO.coefficients[s]=sumO.coefficients[s]+alpha_aos[alpha_ao_index];
				  }
			  }
			  //apply gamma
			  for (int s=0;s<vecStateEnum.size();s++) {
				  sumO.coefficients[s]=sumO.coefficients[s]*gamma;
			  }
			  //add reward
			  for (int s=0;s<vecStateEnum.size();s++) {
				  pair<int,string> rewardInput{s,a};
				  sumO.coefficients[s]=sumO.coefficients[s]*reward[rewardInput];
			  }
			  v_alpha_ab.push_back(sumO); //add the new vector to the list
		  }
	}
	//get backup value
	int maxIndex=0;
	double maxValue=0;
	//computer argMax alpha_ab  of b*alpha_ab
	for (AlphaVector alpha_ab:v_alpha_ab) {
		double sum=0;
		for (int s=0; s<vecStateEnum.size();s++) {
			sum=sum+b[s]*alpha_ab.coefficients[s];
		}
		if (sum>maxValue) {
			maxValue=sum;
			maxIndex=maxValue;
		}
	}
	return v_alpha_ab[maxIndex];
	}


vector<AlphaVector> Mdp::improveValue(vector<AlphaVector> vAlpha,vector<map<int,double>> bSet) {

	std::map<std::tuple<int,string,int,int>,double> alpha_aos;
	for (int i=0; i<vAlpha.size();i++) {
		alpha_aos=calculateAlphaAOS(i,vAlpha[i],alpha_aos);

	}
	bool updated=true;
	while (updated) {
		updated=false;
		for (auto b:bSet) {
			//calculate alpha_ao to use in the backup function
			AlphaVector alpha=backup(b,vAlpha,alpha_aos);
			auto it=std::find(vAlpha.begin(),vAlpha.end(),alpha);
			if (it==vAlpha.end()) {
				vAlpha.push_back(alpha);
				int i=vAlpha.size()-1;
				//get new alpha_aos
				alpha_aos=calculateAlphaAOS(i,vAlpha[i],alpha_aos);
				updated=true;
			}
		}
	}
	return vAlpha;
}

vector<map<int,double>> Mdp::expandBeliefSet(vector<map<int,double>> beliefSet) {
	vector<map<int,double>> newBeliefSet=beliefSet;
	for (auto b:beliefSet) {
		vector<map<int,double>> updatedBeliefSet;
		for (string a:actions) {
			for (int o=0; o<vecObsEnum.size();o++) {
				double prob=1;
				for (int s=0; s<vecStateEnum.size();s++) {
					auto obsIndex=std::make_tuple(o,s,a);
					prob=prob*observationFunction[obsIndex];
				}
				if (prob>0) {
					VariableSet observations=vecObsEnum[o];
				    map<int,double> newB=pomdpUpdateBelief(b,a, std::map<string, string>(), observations.set);
				    updatedBeliefSet.push_back(newB);
				}
			}
		}
		int maxIndex;
		double maxL=-1000;
		for (int i=0; i<updatedBeliefSet.size();i++) {
			double l=0;
			for (int s=0;s<vecStateEnum.size();s++) {
				l=l+b[s]-updatedBeliefSet[i][s];
			}
			if (l>maxL) {
				maxL=l;
				maxIndex=i;
			}
		}
		newBeliefSet.push_back(updatedBeliefSet[maxIndex]);
	}
	return newBeliefSet;
}


void Mdp::pomdpValueIteration(bool rewrite) {
	std::vector<std::map<int,double>> bSet;
	bSet.push_back(this->belief);
	double maxDiff=0;
	double epsilon=0.1;


	AlphaVector a_0;
	a_0.action=this->actions[0];
	a_0.coefficients.insert(a_0.coefficients.begin(), vecStateEnum.size(),5);
	vector<AlphaVector> vAlpha,vAlphaNew;
	vAlpha.push_back(a_0);
	int i=0;
	int max=12;
    do {
    	vAlphaNew=improveValue(vAlpha,bSet);
    	bSet=expandBeliefSet(bSet);

    	double maxDiff=0;


        } while (i<max);

}



//UTILITY FUNCTIONS

void Mdp::printTransitionFunction() {
    cout << "Transition Function\n";
    //    for (int i = 0; i < vecStateEnum.size(); i++) {
    for (string action : actions) {
        pair<int, string> tInput{528, action};
        std::map<int, double> tOutput = transition[tInput];

        VariableSet vs = vecStateEnum[528];
        for (auto s : vs.set) {
            cout << s.first << " " << s.second << " ";
        }
        cout << "\n";
        cout << action << "\n";
        for (auto out : tOutput) {
            VariableSet vo = vecStateEnum[out.first];
            for (auto s : vo.set) {
                cout << s.first << " " << s.second << " ";
            }
            cout << "\n";
            cout << out.first << "\n";
        }
        cout << "\n";
    }
    //    }
}

void Mdp::printRewardFunction() {
    cout << "Reward Function\n";
    for (auto el : reward) {
        if (el.second > 0) {
            VariableSet state = vecStateEnum[el.first.first];
            for (auto s : state.set) {

                cout << s.first << " " << s.second << "\n";
            }
            cout << el.first.second << " " << el.second << "\n";
        }
    }
}

void Mdp::printStates() {
    cout << "\n\nPrint Variables\n";
    for (std::map<string, std::vector<string> >::iterator i = varValues.begin(); i != varValues.end(); i++) {
        cout << i->first << "\n";
        cout << "Values= ";
        for (int j = 0; j < i->second.size(); j++) {
            cout << i->second[j] << " ";
        }
        cout << "\n";
    }
}

void Mdp::printBelief() {
    int k = 0;
    for (auto b : belief) {
        VariableSet vs = vecStateEnum[b.first];
        k++;
        cout << "belief " << k << "\n";
        cout << "probability " << b.second << "\n";
        for (auto j : vs.set) {
            cout << j.first << " " << j.second << " ";
        }
        cout << "\n";
    }
}

void Mdp::createMdp(string name, bool rewrite) {
    this->name = name;

    string fileName = name + ".pomdp";

    std::vector<std::vector<int>> enumInput;
    for (string variable : variables) {
        std::vector<int> valValues;
        for (int i = 0; i < varValues[variable].size(); i++) {
            valValues.push_back(i);
        }
        enumInput.push_back(valValues);
    }

    NestedLoop loop(enumInput);
    std::vector<std::vector<int>> enumOutput = loop.buildMatrix();
    for (int i = 0; i < enumOutput.size(); i++) {
        VariableSet v;
        for (int j = 0; j < enumOutput[i].size(); j++) {
            string name = variables[j];
            std::vector<string> values = varValues[name];
            v.set[name] = values[enumOutput[i][j]];
        }

        mapStateEnum[v] = i;
        //        if (v.set["cup_capacity"] == "0" && v.set["cup_containsWater"] == "0" && v.set["cup_isAt"] == "table" && v.set["cup_isHot"] == "0"
        //                && v.set["glass_capacity"] == "1" && v.set["glass_containsWater"] == "1" && v.set["glass_isAt"] == "table" && v.set["glass_isHot"] == "0"
        //                && v.set["human_isAt"] == "table" && v.set["waterBottle_capacity"] == "0" && v.set["waterBottle_isAt"] == "human") {
        //            cout << "ah" << "\n";
        //            int k = mapStateEnum[v];
        //            cout << k << "\n";
        //        }
        vecStateEnum.push_back(v);
    }

    ifstream inputFile(fileName);
    if (inputFile.good() && !rewrite) {
        for (int i = 0; i < vecStateEnum.size(); i++) {
            for (string action : actions) {
                pair<int, string> transitionInput{i, action};
                std::map<int, double> transitionOutput;


                string line;
                getline(inputFile, line);
                int i1, i2;
                i1 = 0;

                pair<int, string> bTransitionInput{i, action};

                i1 = i2 + 1;
                while ((i2 = line.find(" ", i1)) != string::npos) {
                    string s = line.substr(i1, i2 - i1);
                    i1 = i2 + 1;
                    i2 = line.find(" ", i1);
                    string b = line.substr(i1, i2 - 1);
                    i1 = i2 + 1;
                    transitionOutput[stoi(s)] = stod(b);

                    std::vector<int> previousBeliefs = predecessors[bTransitionInput];
                    previousBeliefs.push_back(stoi(s));
                    predecessors[bTransitionInput] = previousBeliefs;

                }
                //                string b = line.substr(i1);
                //                transitionOutput[stoi(s)] = stod(b);

                transition[transitionInput] = transitionOutput;

                getline(inputFile, line);
                pair<int, string> rewardInput = {i, action};
                reward[rewardInput] = stoi(line);
            }
        }
        inputFile.close();
    } else {
        ofstream file(fileName);


        cout << "Starting Enumeration\n";
        for (int i = 0; i < vecStateEnum.size(); i++) {
            if (isGoalState(vecStateEnum[i])) {
                goalStates.push_back(i);
            }

            for (string action : actions) {
                std::map<VariableSet, double>futureBeliefs = transitionFunction(vecStateEnum[i], action);
                std::map<int, double> transitionOutput;

                pair<int, string> transitionInput{i, action};
                for (auto belief : futureBeliefs) {

                    //                    if (i == 528 && action == "human_fill_glass_waterBottle") {
                    //                        for (auto s : belief.first.set) {
                    //                            cout << s.first << " " << s.second << "\n";
                    //                        }
                    //                        
                    //                        VariableSet tv=vecStateEnum[578];
                    //                        for (auto s:tv.set) {
                    //                            cout<<s.first<<" "<<s.second<<"\n";
                    //                        }
                    //                    }


                    int s = mapStateEnum[belief.first];

                    transitionOutput[s] = belief.second;

                    file << s << " " << belief.second << " ";

                    pair<int, string> bTransitionInput{s, action};
                    std::vector<int> previousBeliefs = predecessors[bTransitionInput];
                    previousBeliefs.push_back(i);
                    predecessors[bTransitionInput] = previousBeliefs;
                }
                file << "\n";
                transition[transitionInput] = transitionOutput;



                pair<int, string> rewardInput{i, action};
                reward[rewardInput] = rewardFunction(vecStateEnum[i], action);
                file << reward[rewardInput] << "\n";
            }
        }
        file.close();
        //        cout << predecessors.size() << "\n";
    }
}

string Mdp::chooseAction() {
    string bestAction;
    double maxValue = -1;
    std::map<string, int> actionValues;
    for (auto aBelief : belief) {
        int s = aBelief.first;
        double prob = aBelief.second;
        for (string action : actions) {
            pair<int, string> qInput{s, action};
            actionValues[action] = actionValues[action] + qValue[qInput] * prob;
        }
    }
    for (auto av : actionValues) {
        if (av.second > maxValue) {
            maxValue = av.second;
            bestAction = av.first;
        }
    }

    return bestAction;

}

void Mdp::setInitialState(std::vector<pair<VariableSet, double>> initialBelief) {
    for (auto aBelief : initialBelief) {
        int s = mapStateEnum[aBelief.first];
        belief[s] = aBelief.second;
    }
}

double Mdp::getActionValue(string action) {
    int actionValue = 0;
    int i = 0;
    for (auto b : belief) {
        i++;
        actionValue = actionValue + qValue[pair<int, string>(b.first, action)] * b.second;
    }
    return actionValue / i;
}

void Mdp::printActualQValues() {
    std::map<string, int> actionValues;
    for (auto aBelief : belief) {
        int s = aBelief.first;
        double prob = aBelief.second;
        for (string action : actions) {
            pair<int, string> qInput{s, action};
            actionValues[action] = actionValues[action] + qValue[qInput] * prob;
        }
    }
    cout << "Q-Values\n";
    for (auto av : actionValues) {
        cout << av.first << " " << av.second << "\n";
    }

}

std::map<string, double> Mdp::getNormQValue() {
    std::map<string, double> result, result2;

    for (auto aBelief : belief) {
        int s = aBelief.first;
        double prob = aBelief.second;
        for (string action : actions) {
            pair<int, string> qInput{s, action};
            result[action] = result[action] + qValue[qInput] * prob;
        }
    }
    double sum = 0;
    for (auto r : result) {
        sum = sum + r.second;
    }
    for (auto r : result) {
        result2[r.first] = r.second / sum;
    }
    return result2;
}

double Mdp::getQValue(string action) {
    double result;
    for (auto aBelief : belief) {
        int s = aBelief.first;
        double prob = aBelief.second;
        pair<int, string> qInput{s, action};
        result = result + qValue[qInput] * prob;
    }
}
