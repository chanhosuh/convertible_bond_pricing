/* 
 * File:   CIRbond.cpp
 * Author: suh
 * 
 * Created on May 4, 2013, 2:01 AM
 */

#include "CIRbond.h"
#include <cmath>

CIRbond::CIRbond(double mean, double speed, double maturity, double sigma)
				: mean(mean), speed(speed), maturity(maturity), sigma(sigma)
{
	gamma = 0.5*sqrt(speed*speed + 2*sigma*sigma);
}

double CIRbond::C(double t)
{
	double a = speed;
	double b = mean;
	double T = maturity;

	return sinh(gamma*(T-t))/(gamma*cosh(gamma*(T-t)) + 0.5*a*sinh(gamma*(T-t)));
}

double CIRbond::A(double t)
{
	double a = speed;
	double b = mean;
	double T = maturity;

	return -2*a*b/(sigma*sigma)*log(gamma*exp(0.5*a*(T-t))/(gamma*cosh(gamma*(T-t)) + 0.5*a*sinh(gamma*(T-t))));
}

double CIRbond::price(double r, double t)
{
	return exp(-r*C(t) - A(t));
}