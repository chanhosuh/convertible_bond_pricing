/* 
 * File:   CIRbond.h
 * Author: suh
 *
 * Created on May 4, 2013, 2:01 AM
 */

#ifndef _CIRBOND_H
#define	_CIRBOND_H

class CIRbond {
public:
    CIRbond(double mean, double speed, double maturity, double sigma);
    virtual ~CIRbond() {}

	virtual double price(double r, double t);

private:
	double mean;
	double speed;
	double maturity;
	double sigma;
	double gamma;

	double A(double t);
	double C(double t);

};

#endif	/* _CIRBOND_H */

