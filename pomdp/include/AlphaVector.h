/*
 * AlphaVector.h
 *
 *  Created on: May 9, 2015
 *      Author: mfiore
 */

#ifndef SOURCE_DIRECTORY__UTILITIES_POMDP_ALPHAVECTOR_H_
#define SOURCE_DIRECTORY__UTILITIES_POMDP_ALPHAVECTOR_H_

#include <string>
#include <vector>

using namespace std;
class AlphaVector {
public:
	AlphaVector();
	AlphaVector(const AlphaVector & other);

	inline bool operator==(const AlphaVector& rhs) {
		if (this->coefficients.size() != rhs.coefficients.size()) return false;
		for (int i=0; i<this->coefficients.size(); i++) {
			if (this->coefficients[i]!=rhs.coefficients[i]) return false;
		}
		return this->action==rhs.action;
	}
	virtual ~AlphaVector();

	vector<double> coefficients;
	string action;
};

#endif /* SOURCE_DIRECTORY__UTILITIES_POMDP_ALPHAVECTOR_H_ */
