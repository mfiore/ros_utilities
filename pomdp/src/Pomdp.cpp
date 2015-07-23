/*
 * Pomdp.cpp
 *
 *  Created on: May 9, 2015
 *      Author: mfiore
 */

#include "Pomdp.h"

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <stack>
#include <utility>
#include <algorithm>
#include <c++/4.6/bits/stl_queue.h>

Pomdp::Pomdp() {
}

Pomdp::Pomdp(const Pomdp& orig) {
}

Pomdp::~Pomdp() {
}



/* This function filters the current belief with the observed states, erasing the belief that don't respect them.
 * obsStates: the observed states and their values
 */
map<int,double> Pomdp::filterBelief(map<int,double> belief,std::map<string, string> obsStates) {
    std::map<int, double> newBelief;
    bool checkObservedVariables = true;

    for (auto beliefState : belief) {
        checkObservedVariables = true;

        //we check if all the observed variables in the belief are equal the observed variables in input
        VariableSet v = vecStateEnum_[beliefState.first];
        for (auto j : obsStates) {
            if (v.set.at(j.first) != j.second) {
                checkObservedVariables = false;
                break;
            }
        }
        //if so we keep the belief
        if (checkObservedVariables) {
            newBelief[beliefState.first] = beliefState.second;

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


map<int,double> Pomdp::updateBelief(map<int,double> belief,string action, std::map<string, string> obsStates, int o) {
	  std::map<int, double> finalBelief;

	    //for each current belief
	    for (auto beliefState: belief) {
	        int s = beliefState.first;

	        pair<int, string> transitionInput{s, action};
	        std::map<int, double> futureBelief = transitionMap_[transitionInput];

	        for (auto futureBeliefState: futureBelief) {
	            finalBelief[futureBeliefState.first] = finalBelief[futureBeliefState.first] + futureBeliefState.second;
	        }

	        //normalize
	        double sum = 0;
	        for (auto beliefState : finalBelief) {
	            sum = sum + beliefState.second;
	        }
	        for (auto beliefState: finalBelief) {
	            beliefState.second = beliefState.second / sum;
	        }
	    }
	    //apply observation function
	    for (auto beliefState:finalBelief) {
	    	std::tuple<int,int,string> observationInput{o,beliefState.first,action};
	    	beliefState.second=beliefState.second*observationMap_[observationInput];
	    }
	    //normalize
	     double sum = 0;
	     for (auto beliefState : finalBelief) {
	         sum = sum + beliefState.second;
	     }
	     for (auto beliefState: finalBelief) {
	         beliefState.second = beliefState.second / sum;
	     }
	  	finalBelief=filterBelief(finalBelief,obsStates);

	    return finalBelief;
}


/*
 This function updates the current system belief
 *
 * action: current action
 * obsStates: set of observed variables
 * observations: set of observations
 */
void Pomdp::updateCurrentBelief(string action, std::map<string, string> obsStates, std::map<string, string> observations) {
	map<int,double> newBelief;

	VariableSet vs;
	vs.set=observations;
	int o=mapObservationEnum_[vs];

	newBelief=updateBelief(this->currentBelief_,action,obsStates,o);
    this->currentBelief_=newBelief;
}



std::map<std::tuple<int,string,int,int>,double> Pomdp::calculateAlphaAOS(int i, AlphaVector alpha,
		std::map<std::tuple<int,string,int,int>,double> currentSet) {
	for (string a:actions_) {
		for (int o=0; o<vecObsEnum_.size(); o++) {
			for (int s=0; s<vecStateEnum_.size(); s++) {
				double sum=0;
				for (int sp=0;s<vecStateEnum_.size();s++) {
					std::tuple<int,int,string> observationInput=std::make_tuple(o,sp,a);
					pair<int,string> transitionInput{s,a};
					sum=sum+alpha.coefficients[sp]*observationMap_[observationInput]*transitionMap_[transitionInput][sp];
				}
				auto alphaIndex=std::make_tuple(i,a,o,s);
				currentSet[alphaIndex]=sum;
			}
		}
	}
	return currentSet;
}

AlphaVector Pomdp::backup(map<int,double> b,vector<AlphaVector> vAlpha,
		std::map<std::tuple<int,string,int,int>,double> alpha_aos) {

	double gamma=0.3;
	//compute alpha_ab
	vector<AlphaVector> v_alpha_ab;
	//there is a different alpha_ab for each alpha and for each a
	for (string a:actions_) {
		  for (AlphaVector alpha:vAlpha) {
			  AlphaVector sumO;
			  sumO.action=a;
			  for (int o=0; o<vecObsEnum_.size();o++) {
				  //find argmax a_o| a_o*b
				  double max_alpha_ao;
				  int max_alpha_index=0;
				  for (int alphaIndex=0; alphaIndex<vAlpha.size();alphaIndex++) {
					  //get this alpha_ao
					  double sum=0;
					  for (int s=0; s<vecStateEnum_.size();s++) {
						  auto alpha_ao_index=std::make_tuple(alphaIndex,a,o,s);
						  sum=b[s]*alpha_aos[alpha_ao_index];
					  }
					  if (sum>max_alpha_ao) {
						  max_alpha_ao=sum;
						  max_alpha_index=alphaIndex;
					  }
				  }
				  //sum argMax to current max sum
				  for (int s=0; s<vecStateEnum_.size();s++) {
					  auto alpha_ao_index=std::make_tuple(max_alpha_index,a,o,s);
					  sumO.coefficients[s]=sumO.coefficients[s]+alpha_aos[alpha_ao_index];
				  }
			  }
			  //apply gamma
			  for (int s=0;s<vecStateEnum_.size();s++) {
				  sumO.coefficients[s]=sumO.coefficients[s]*gamma;
			  }
			  //add reward
			  for (int s=0;s<vecStateEnum_.size();s++) {
				  pair<int,string> rewardInput{s,a};
				  sumO.coefficients[s]=sumO.coefficients[s]*rewardMap_[rewardInput];
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
		for (int s=0; s<vecStateEnum_.size();s++) {
			sum=sum+b[s]*alpha_ab.coefficients[s];
		}
		if (sum>maxValue) {
			maxValue=sum;
			maxIndex=maxValue;
		}
	}
	return v_alpha_ab[maxIndex];
	}


vector<AlphaVector> Pomdp::improveValue(vector<AlphaVector> vAlpha,vector<map<int,double>> bSet) {

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

vector<map<int,double>> Pomdp::expandBeliefSet(vector<map<int,double>> beliefSet) {
	vector<map<int,double>> newBeliefSet=beliefSet;
	for (auto b:beliefSet) {
		vector<map<int,double>> updatedBeliefSet;
		for (string a:actions_) {
			for (int o=0; o<vecObsEnum_.size();o++) {
				double prob=1;
				for (int s=0; s<vecStateEnum_.size();s++) {
					auto obsIndex=std::make_tuple(o,s,a);
					prob=prob*observationMap_[obsIndex];
				}
				if (prob>0) {
					map<int,double> newB=updateBelief(b,a, std::map<string, string>(), o);
				    updatedBeliefSet.push_back(newB);
				}
			}
		}
		int maxIndex;
		double maxL=-1000;
		for (int i=0; i<updatedBeliefSet.size();i++) {
			double l=0;
			for (int s=0;s<vecStateEnum_.size();s++) {
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


void Pomdp::valueIteration(bool rewrite) {
	std::vector<std::map<int,double>> bSet;
	bSet.push_back(this->currentBelief_);
	double maxDiff=0;
	double epsilon=0.1;


	AlphaVector a_0;
	a_0.action=this->actions_[0];
	a_0.coefficients.insert(a_0.coefficients.begin(), vecStateEnum_.size(),5);
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

void Pomdp::printTransitionFunction() {
    cout << "Transition Function\n";
    //    for (int i = 0; i < vecStateEnum.size(); i++) {
    for (string action : actions_) {
        pair<int, string> tInput{528, action};
        std::map<int, double> tOutput = transitionMap_[tInput];

        cout << action << "\n";
        for (auto out : tOutput) {
            VariableSet vo = vecStateEnum_[out.first];
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

void Pomdp::printRewardFunction() {
    cout << "Reward Function\n";
    for (auto el : rewardMap_) {
        if (el.second > 0) {
            VariableSet state = vecStateEnum_[el.first.first];
            for (auto s : state.set) {

                cout << s.first << " " << s.second << "\n";
            }
            cout << el.first.second << " " << el.second << "\n";
        }
    }
}

void Pomdp::printStates() {
    cout << "\n\nPrint Variables\n";
    for (std::map<string, std::vector<string> >::iterator i = varValues_.begin(); i != varValues_.end(); i++) {
        cout << i->first << "\n";
        cout << "Values= ";
        for (int j = 0; j < i->second.size(); j++) {
            cout << i->second[j] << " ";
        }
        cout << "\n";
    }
}

void Pomdp::printBelief() {
    int k = 0;
    for (auto b : currentBelief_) {
        VariableSet vs = vecStateEnum_[b.first];
        k++;
        cout << "belief " << k << "\n";
        cout << "probability " << b.second << "\n";
        for (auto j : vs.set) {
            cout << j.first << " " << j.second << " ";
        }
        cout << "\n";
    }
}

void Pomdp::createPomdp(string name, bool rewrite) {
    this->name_ = name;

    string fileName = name + ".pomdp";

    std::vector<std::vector<int>> enumInput;
    for (string variable : variables_) {
        std::vector<int> valValues;
        for (int i = 0; i < varValues_[variable].size(); i++) {
            valValues.push_back(i);
        }
        enumInput.push_back(valValues);
    }

    NestedLoop loop(enumInput);
    std::vector<std::vector<int>> enumOutput = loop.buildMatrix();
    for (int i = 0; i < enumOutput.size(); i++) {
        VariableSet v;
        for (int j = 0; j < enumOutput[i].size(); j++) {
            string name = variables_[j];
            std::vector<string> values = varValues_[name];
            v.set[name] = values[enumOutput[i][j]];
        }

        mapStateEnum_[v] = i;
        vecStateEnum_.push_back(v);
    }

    enumInput.clear();
    for (string observation:observationsVariables_) {
    	  std::vector<int> valValues;
    	        for (int i = 0; i < observationValues_[observation].size(); i++) {
    	            valValues.push_back(i);
    	        }
    	        enumInput.push_back(valValues);
    }
    NestedLoop loop2(enumInput);
    enumOutput = loop2.buildMatrix();
    for (int i = 0; i < enumOutput.size(); i++) {
        VariableSet v;
        for (int j = 0; j < enumOutput[i].size(); j++) {
            string name = variables_[j];
            std::vector<string> values = observationValues_[name];
            v.set[name] = values[enumOutput[i][j]];
        }

        mapObservationEnum_[v] = i;
        vecObsEnum_.push_back(v);
    }


    ifstream inputFile(fileName);
    if (inputFile.good() && !rewrite) {
        for (int i = 0; i < vecStateEnum_.size(); i++) {
            for (string action : actions_) {
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

                }
                //                string b = line.substr(i1);
                //                transitionOutput[stoi(s)] = stod(b);

                transitionMap_[transitionInput] = transitionOutput;

                getline(inputFile, line);
                pair<int, string> rewardInput = {i, action};
                rewardMap_[rewardInput] = stoi(line);
            }
        }
        inputFile.close();
    } else {
        ofstream file(fileName);


        cout << "Starting Enumeration\n";
        for (int i = 0; i < vecStateEnum_.size(); i++) {

            for (string action : actions_) {
                std::map<VariableSet, double>futureBelief = transitionFunction(vecStateEnum_[i], action);
                std::map<int, double> transitionOutput;

                pair<int, string> transitionInput{i, action};
                for (auto beliefState : futureBelief) {

                    int s = mapStateEnum_[beliefState.first];

                    transitionOutput[s] = beliefState.second;

                    file << s << " " << beliefState.second << " ";

                }
                file << "\n";
                transitionMap_[transitionInput] = transitionOutput;

                for (int o=0; o<vecObsEnum_.size(); o++) {
                	VariableSet observationState=vecObsEnum_[i];
                	auto observationInput=std::make_tuple(o,i,action);
                	observationMap_[observationInput]=observationFunction(observationState,vecStateEnum_[i],action);

                	file<< o << " "  << i << " " << action << observationMap_[observationInput]<<"\n";
                }
                file<<"\n";

                pair<int, string> rewardInput{i, action};
                rewardMap_[rewardInput] = rewardFunction(vecStateEnum_[i], action);
                file << rewardMap_[rewardInput] << "\n";
            }
        }
        file.close();
        //        cout << predecessors.size() << "\n";
    }
}

string Pomdp::chooseAction() {
	double max=0;
	int max_index=0;
	for (int i=0; i<alphaVectors_.size();i++) {
		double sum=0;
		for (int s=0; s<vecStateEnum_.size();s++) {
			sum=sum+currentBelief_[s]*alphaVectors_[i].coefficients[s];
		}
		if (sum>max) {
			max=sum;
			max_index=i;
		}
	}
    return alphaVectors_[max_index].action;

}

void Pomdp::setInitialState(std::vector<pair<VariableSet, double>> initialBelief) {
    for (auto beliefState : initialBelief) {
        int s = mapStateEnum_[beliefState.first];
        currentBelief_[s] = beliefState.second;
    }
}


