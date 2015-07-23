/* 
 * File:   PriorityQueue.cpp
 * Author: mfiore
 * 
 * Created on January 28, 2015, 10:17 AM
 */

#include "PriorityQueue.h"

PriorityQueue::PriorityQueue() {
}

PriorityQueue::PriorityQueue(const PriorityQueue& orig) {
}

PriorityQueue::~PriorityQueue() {
}

bool cmp(pair<int, double>* left, pair<int, double>* right) {
    return left->second < right->second;
}

bool PriorityQueue::isEmpty() {
    return queue.empty();
}

void PriorityQueue::pushElement(pair<int, double> element) {
    if (indexs.find(element.first) != indexs.end()) {
        editElement(element.first, element.second);
    } else {
        pair<int, double> * newElement = new pair<int, double> (element.first, element.second);

        if (queue.empty()) {
            queue.push_back(newElement);
            std::make_heap(queue.begin(), queue.end(), cmp);

            indexs[element.first] = newElement;
        } else {
            queue.push_back(newElement);
            std::push_heap(queue.begin(), queue.end(), cmp);
            indexs[element.first] = newElement;

        }
    }
}

void PriorityQueue::editElement(int index, double priority) {
    pair<int, double> * element = indexs[index];
    element->second = priority;
    std::make_heap(queue.begin(), queue.end(), cmp);
}

pair<int, double> PriorityQueue::pop() {
    pair<int, double>* element = queue.back();
    pair<int, double> returnElement{element->first, element->second};
    queue.pop_back();
    std::pop_heap(queue.begin(), queue.end(), cmp);
    indexs.erase(element->first);
    delete element;
    return returnElement;
}

double PriorityQueue::getPriority(int i) {
    if (indexs.find(i) != indexs.end()) {
        return indexs[i]->second;
    } else return 0;
}

int PriorityQueue::size() {
    return queue.size();
}

void PriorityQueue::print() {
    for (auto element:queue) {
        cout<<element->first<<" "<<element->second<<"\n";
    }
}


