/* 
 * File:   couponBond.cpp
 * Author: suh
 * 
 * Created on May 4, 2013, 2:20 AM
 */

#include "CouponBond.h"

CouponBond::CouponBond(double par, double couponRate,
					   double mean, double speed, double maturity, double sigma)
	: CIRbond(mean, speed, maturity, sigma), par(par), couponRate(couponRate)
{
}

double CouponBond::price(double r, double t)
{
	return par*CIRbond::price(r, t);
}

