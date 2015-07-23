/*
 * AlphaVector.cpp
 *
 *  Created on: May 9, 2015
 *      Author: mfiore
 */

#include "AlphaVector.h"

AlphaVector::AlphaVector() {
	// TODO Auto-generated constructor stub

}

AlphaVector::~AlphaVector() {
	// TODO Auto-generated destructor stub
}

AlphaVector::AlphaVector(const AlphaVector & other) {
	this->coefficients=other.coefficients;
	this->action=other.action;
}
