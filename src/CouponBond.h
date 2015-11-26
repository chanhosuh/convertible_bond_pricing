/* 
 * File:   couponBond.h
 * Author: suh
 *
 * Created on May 4, 2013, 2:20 AM
 */

#ifndef _COUPONBOND_H
#define	_COUPONBOND_H

#include "CIRbond.h"

class CouponBond : public CIRbond
{
public:
    CouponBond(double par, double couponRate,
			double mean, double speed, double maturity, double sigma);
    virtual ~CouponBond() {}

	virtual double price(double r, double t);

private:
	double par;
	double couponRate;
};

#endif	/* _COUPONBOND_H */

