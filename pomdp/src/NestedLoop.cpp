/* 
 * File:   NestedLoop.cpp
 * Author: mfiore
 * 
 * Created on December 16, 2014, 10:46 AM
 */

#include <queue>
#include <iostream>
#include "NestedLoop.h"

NestedLoop::NestedLoop(vector<vector<int> > variables) {
    this->values = variables;
}

NestedLoop::NestedLoop(const NestedLoop& orig) {
}

NestedLoop::~NestedLoop() {
}

/*
 In this function we want to build a matrix of possible assignment of the variables
 * ex. 0 0 0
 *     0 0 1
 *     0 1 0 ... etc.
 * 
 * and return it in matrix form
 */
vector<vector<int> > NestedLoop::buildMatrix() {
    typedef vector<int> assignment; //just a convenience typedef

    vector<assignment> result; //result will contain the final matrix

    int depthLevel = values.size() - 1; //the number of nested loops (counting from zero)

    //we represent the recursion as a queue which contains the current variable assignment (ex. 0 0) and the current depth level
    std::queue<pair<assignment, int> > q;

    std::cout << values.size() << "\n";
    //we start assigning the external loop values
    for (int i = 0; i < values[0].size(); i++) {
        assignment a;
        a.push_back(values[0][i]);
        pair<assignment, int> p;
        p.first = a;
        p.second = 0;
        q.push(p);
    }
    int i = 0;
    while (!q.empty()) {

        pair<assignment, int> element = q.front();
        i = element.second;
        assignment currentAssignment = element.first;
        q.pop();
     
        /*when i is different then depthLevel we still haven't reached the innermost loop so we just expand the next inner variable 
        ex: if our current element is 0, we could expand this assignment to 0 0 and 0 1 , putting them in the queue. */
        if (i != depthLevel) {
            for (int j = 0; j < values[i + 1].size(); j++) {
                assignment a = currentAssignment;
                a.push_back(values[i + 1][j]);
                pair<assignment, int> newElement;
                newElement.first = a;
                newElement.second = i + 1;
                q.push(newElement);
            }
        }/* if i==depthLevel we reached the innermost cycle, so we copy the assignment in the returning matrix*/
        else {
            result.push_back(currentAssignment);
        }
    }
    return result;
}

