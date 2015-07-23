/* 
 * File:   NestedLoop.h
 * Author: mfiore
 *
 * Created on December 16, 2014, 10:46 AM
 */

#ifndef NESTEDLOOP_H
#define	NESTEDLOOP_H

#include <vector>
#include <stack>

using namespace std;

class NestedLoop {
public:
    NestedLoop(vector<vector<int> > values);
    NestedLoop(const NestedLoop& orig);
    virtual ~NestedLoop();
    
    vector<vector<int> >  buildMatrix();
private:
    vector<vector<int> > values;
    
};

#endif	/* NESTEDLOOP_H */

