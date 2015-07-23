/* 
 * File:   VariableSet.h
 * Author: mfiore
 *
 * Created on August 19, 2014, 1:04 PM
 */

#ifndef WORLDSTATE_H
#define	WORLDSTATE_H

#include <vector>
#include <map>
#include <iostream>
using namespace std;

class VariableSet {
public:
    VariableSet();
    VariableSet(const VariableSet& orig);
    virtual ~VariableSet();

    std::map<string, string> set;

    inline bool operator<(const VariableSet& b) const {

        std::map<string, string>::const_iterator i;
        i = set.begin();
        while (i != set.end()) {
            if (b.set.find(i->first) != b.set.end()) {
                string v1=b.set.at(i->first);
                string v2=i->second;
                if (v1 == "*" || v2 == "*") {
                    i++;
                } else if (v1 < v2) {
                    return true;
                } else if (v1 > v2) {
                    return false;
                } else {
                    i++;
                }
            } else {
                i++;
            }
        }
        return false;

    }

    inline bool contained(const VariableSet &b) const {
        std::map<string, string>::const_iterator i;
        i = set.begin();
        while (i != set.end()) {
            if (b.set.find(i->first) != b.set.end()) {
                if (i->second != b.set.at(i->first) && i->second != "*" && b.set.at(i->first) != "*") return false;
                i++;
            } else return false;
        }
        return true;
    }

    inline bool operator==(const VariableSet& b) const {
        if (set.size() != b.set.size()) return false;
        std::map<string, string>::const_iterator i;
        i = set.begin();
        while (i != set.end()) {
            if (b.set.find(i->first) != b.set.end()) {
                if (i->second != b.set.at(i->first) && i->second != "*" && b.set.at(i->first) != "*") return false;
                i++;
            } else return false;
        }

    }

    inline void operator=(const VariableSet &b) {
        if (b.set.size() > 0) {
            this->set = b.set;
        }
    }

    void print() {
        for (auto s : set) {
            cout << s.first << " " << s.second << "\n";
        }
    }

private:

};

#endif	

