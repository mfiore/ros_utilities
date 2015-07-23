/* 
 * File:   PriorityQueue.h
 * Author: mfiore
 *
 * Created on January 28, 2015, 10:17 AM
 */

#ifndef PRIORITYQUEUE_H
#define	PRIORITYQUEUE_H

#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;



class PriorityQueue {
public:
    PriorityQueue();
    PriorityQueue(const PriorityQueue& orig);
    virtual ~PriorityQueue();

    bool isEmpty();
    void pushElement(pair<int, double> element);
    void editElement(int index, double priority);
    double getPriority(int index);
    pair<int, double> pop();
    int size();

    void print();

private:
    std::vector<pair<int, double>*> queue;
    std::map<int, pair<int, double>*> indexs;

};

#endif	/* PRIORITYQUEUE_H */

