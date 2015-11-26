/*
 * File:   main.cpp
 * Author: suh
 *
 * Created on May 4, 2013, 9:36 PM
 *
 * Version 3: Contingent conversion is handled as conv bond with no
 * barrier with a knock-out option subtracted.
 * Also, changed from log spot to regular spot.
 */

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "CouponBond.h"

using namespace std;

// intensity for default process;
double lambda(double s)
{
return 1/(s*s);
}

int main()
{
    // constants for CIR interest rate process
    double a = 2.0; // mean reversion speed
    double b = 0.05; // mean interest rate
    double sigma_r = 0.4;

    // equity process
    double q = 0.01; //dividend yield;
    double sigma = 0.3;

    double rho = -0.5; // correlation coef between spot and interest rates

	// pure bond parameters
    double T = 1;  // maturity
    double F = 1000; // par value
    double c = 0.02; // coupon rate

	CouponBond myBond(F, c, b, a, T, sigma_r);

    // convertible bond parameters
    double kappa = 10; // conversion ratio
    double callValue = 110;
	double putValue = 90;
    double contingentHurdle = 140; //knock-in trigger for conversion option
    double callHurdle = 140;
    double R_v = 0.4*F; // value recovered from bond default


	// create grid
    double s_max = 300;
    double r_max = 1;
    int I = 50; // number of equity steps
    int J = 50; // number of interest rate steps
    int K = 252; // number of temporal steps
    double v[I+1][J+1];
    double vb[I+1][J+1];

    double deltaT = T/K;
    double deltaR = r_max/J;

    double deltaS = s_max/I;

	int Ib = I; // reset this below: largest index before exceeding cont hurdle

	/* Set PAYOFF terminal values on grid */

	// payoff for convertible with no barrier
    for (int i = 0; i < I+1; ++i)
    {
        for (int j = 0; j < J+1; ++j)
        {
			v[i][j] = max(kappa*i*deltaS, F);
        }
    }

	// due to the singularity at s=0, t=T
	// we use a cubic spline to interpolate between par and recovery values
	// in the payoff.
	// Since we use vanishing delta at s=0, and the call payoff is just F nearby
	// the derivatives of the spline should be 0 at the two points.
	for (int j=0; j < J+1; ++j)
	{
		v[0][j] = R_v;
		v[1][j] = 0.5*(v[0][j] + v[2][j]); // for a uniform grid, the symmetry fits the cubic spline at the average

	}

	// conversion option with knock-out barrier
	for (int i = 0; i < I+1; ++i)
	{
		if (i*deltaS >= contingentHurdle)
		{
			Ib = i;
			for (int j=0; j < J+1; ++j)
				vb[i][j] = 0;
			break;
		}
		else
		{
			for (int j=0; j < J+1; ++j)
				vb[i][j] = max(kappa*i*deltaS - F, 0.0);
		}
	}


    ofstream errorOutput("errors.txt");

        /* Begin Yanenko splitting */

    // each splitting stage will be solved using an iterative method
    // in particular, PSOR, so we set the tolerance level here
    double tol = 0.0005;
	double maxIter = 100000;
	double noIter = 0;

    double w[I+1][J+1]; // will contain our half-step values
	double wb[I+1][J+1];
    double temp; // temp values so we can compute the error

    for (int k = 1; k < K + 1; ++k)
    {

        // set initial guess for SOR
        for (int i = 0; i < I + 1; ++i)
            for (int j = 0; j < J + 1; ++j)
			{
                w[i][j] = v[i][j]; // start with last time step values
				wb[i][j] = vb[i][j];
			}

        double error = 100; // set error large so we at least enter while loop

        // first half of splitting
		noIter = 0;
        while (error > tol && noIter++ < maxIter) // keep iterating until error is less than tolerance
        {
            error = 0;
            // now loop through using Gauss-Seidel/SOR
            // to compute k-th time step values
// BOUNDARY
// CONDITIONS
            // set "bottom" boundary values:

            // for equity
            for (int j = 0; j < J + 1; ++j)
            {
                // store old values so we can check the error afterwards
                temp = w[0][j];

                // set new values using boundary conditions
                w[0][j] = w[1][j]; // vanishing delta
				//w[0][j] = R_v;
				error = max(error, abs(temp-w[0][j]));

				temp = wb[0][j];
				wb[0][j] = w[1][j];
				error = max(error, abs(temp-wb[0][j]));

            }
            // for interest rate
            for (int i = 1; i < I; ++i)
            {
                // store old values so we can check the error afterwards
                temp = w[i][0];

                // set new values using boundary conditions
                w[i][0] = 2 * w[i][1] - w[i][2]; // vanishing gamma

                // store biggest error so far;
                error = max(error, abs(temp - w[i][0]));

				temp = wb[i][0];
				wb[i][0] = 2*wb[i][1] - wb[i][2];
				error = max(error, abs(temp - wb[i][0]));
            }


            // INTERIOR

            // compute interior values
            for (int i = 1; i < I; ++i)
			{
					if (i*deltaS > callHurdle)
					{
						for (int j=1; j<J; ++j)
						{
							temp = w[i][j];
							w[i][j] = max(kappa*callValue, kappa*i*deltaS);
						    error = max(error, abs(w[i][j] - temp));
						}
					}
					else
					{
						for (int j = 1; j < J; ++j)
						{
							temp = w[i][j];
							w[i][j] = 0.5*(i*i*sigma*sigma + i*(j*deltaR - q + lambda(i*deltaS) - 0.5*sigma*sigma))*w[i+1][j]
									 +0.5*(i*i*sigma*sigma - i*(j*deltaR - q + lambda(i*deltaS) - 0.5*sigma*sigma))*w[i-1][j]
									 +v[i][j]/deltaT
									 +1/8 * rho *sqrt(j*deltaR)*sigma_r*sigma*i/deltaR *(v[i+1][j+1] + v[i-1][j-1] - v[i+1][j-1] - v[i-1][j+1]);
							w[i][j] = w[i][j] / (1/deltaT + i*i*sigma*sigma);
							w[i][j] = max(kappa*i*deltaS, w[i][j]);
							w[i][j] = max(kappa*putValue, w[i][j]);


							// keep track of largest error so far in this iterative step
							error = max(error, abs(w[i][j] - temp));
						}
					}
			}
			for (int i = 1; i < Ib; ++i)
				for (int j=1; j < J; ++j)
				{
					temp = wb[i][j];
					wb[i][j] = 0.5*(i*i*sigma*sigma + i*(j*deltaR - q + lambda(i*deltaS) - 0.5*sigma*sigma))*wb[i+1][j]
						     +0.5*(i*i*sigma*sigma - i*(j*deltaR - q + lambda(i*deltaS) - 0.5*sigma*sigma))*wb[i-1][j]
						     +vb[i][j]/deltaT
						     +1/8 * rho *sqrt(j*deltaR)*sigma_r*sigma*i/deltaR *(vb[i+1][j+1] + vb[i-1][j-1] - vb[i+1][j-1] - vb[i-1][j+1]);
                    wb[i][j] = wb[i][j] / (1/deltaT + i*i*sigma*sigma);
                    wb[i][j] = max( max(kappa*i*deltaS - w[i][j], 0.0), wb[i][j]);

                    // keep track of largest error so far in this iterative step
                    error = max(error, abs(wb[i][j] - temp));
                }

            // BOUNDARY
            // CONDITIONS
            // set "top" boundary values:
            for (int j = 0; j < J + 1; ++j)
            {
                // store old values so we can check the error afterwards
                temp = w[I][j];

                // set new values using boundary conditions
				//w[I][j] = deltaT*kappa + w[I-1][j];
				w[I][j] = 2*w[I-1][j] - w[I-2][j];

                // store biggest error so far;
                error = max(error, abs(temp - w[I][j]));

            }
            for (int i = 1; i < I; ++i)
            {
                // store old values so we can check the error afterwards
                temp = w[i][J];

                // set new values using boundary conditions
                w[i][J] = w[i][J - 1]; // vanishing rho

                // store biggest error so far;
                error = max(error, abs(temp - w[i][J]));

				temp = wb[i][J];
				wb[i][J] = wb[i][J-1];

				error = max(error, abs(temp - wb[i][J]));
            }

            errorOutput << error << endl;

        }
        // while post-condition: w contains approximate solution
        // which was close enough to last approximation in the iteration.
        // w is the solution at the time step "(k-1) + 1/2)".

        // now that we have half-step values w[i][j],
        // begin second half of split solving

        // set initial guess for SOR
        for (int i = 0; i < I + 1; ++i)
            for (int j = 0; j < J + 1; ++j)
			{
                v[i][j] = w[i][j]; // start with last half-step values
				vb[i][j] = wb[i][j];
			}

        error = 100; // reset error so that we at least go through while loop once

		noIter = 0;
        while (error > tol && noIter++ < maxIter)
        {
            error = 0;

            // begin SOR

            // BOUNDARY
            // CONDITIONS
            // for equity
            for (int j = 0; j < J + 1; ++j)
            {
                // store old values so we can check the error afterwards
                temp = v[0][j];
                //v[k][0][j] = R_v;
                v[0][j] = v[1][j]; // vanishing delta
				//v[0][j] = R_v;
				error = max(error, abs(temp-v[0][j]));


                temp = vb[0][j];
                vb[0][j] = vb[1][j];
				error = max(error, abs(temp-vb[0][j]));

            }
            //for interest rate
            for (int i = 1; i < I; ++i)
            {
                // store old values so we can check the error afterwards
                temp = v[i][0];

                // set new values using boundary conditions
                v[i][0] = 2 * v[i][1] - v[i][2]; // vanishing gamma

                // store biggest error so far;
                error = max(error, abs(temp - v[i][0]));

                temp = vb[i][0];
                vb[i][0] = 2 * vb[i][1] - vb[i][2]; // vanishing gamma

                // store biggest error so far;
                error = max(error, abs(temp - vb[i][0]));
            }



            // INTERIOR
            // compute interior values
            for (int i = 1; i < I; ++i)
			{
				if (i*deltaS > callHurdle)
					for (int j=1; j < J; ++j)
					{
						temp = v[i][j];
						v[i][j] = max(kappa*callValue, kappa*i*deltaS);
						error = max(error, abs(v[i][j] - temp));
					}
				else
						for (int j = 1; j < J; ++j)
						{
							temp = v[i][j];
							v[i][j] = 0.5/deltaR*(j*sigma_r*sigma_r + a*(b-j*deltaR))*v[i][j+1]
									 +0.5/deltaR*(j*sigma_r*sigma_r - a*(b-j*deltaR))*v[i][j-1]
									 +w[i][j]/deltaT
									 +1/8*rho*sqrt(j*deltaR)*sigma_r*sigma*i/deltaR*(w[i+1][j+1]+w[i-1][j-1]-w[i+1][j-1]-w[i-1][j+1])
									 +lambda(i*deltaS)*R_v;
							v[i][j] = v[i][j] / (1/deltaT + j/deltaR*sigma_r*sigma_r + j*deltaR + lambda(i*deltaS));
							v[i][j] = max(kappa*i*deltaS, v[i][j]);
							v[i][j] = max(kappa*putValue, v[i][j]);
							error = max(error, abs(v[i][j] - temp));
						}
			}
			for (int i=1; i < Ib; ++i)
				for (int j=1; j < J; ++j)
				{
					temp = vb[i][j];
					vb[i][j] = 0.5/deltaR*(j*sigma_r*sigma_r + a*(b-j*deltaR))*vb[i][j+1]
							 +0.5/deltaR*(j*sigma_r*sigma_r - a*(b-j*deltaR))*vb[i][j-1]
						     +wb[i][j]/deltaT
						     +1/8*rho*sqrt(j*deltaR)*sigma_r*sigma*i/deltaR*(wb[i+1][j+1]+wb[i-1][j-1]-wb[i+1][j-1]-wb[i-1][j+1]);
                    vb[i][j] = vb[i][j] / (1/deltaT + j/deltaR*sigma_r*sigma_r + j*deltaR + lambda(i*deltaS));
                    vb[i][j] = max( max(kappa*i*deltaS - v[i][j], 0.0), vb[i][j]);

                    error = max(error, abs(vb[i][j] - temp));
                }

            // BOUNDARY
            // CONDITIONS
            // set boundary values:
            // for equity
            for (int j = 0; j < J + 1; ++j)
            {
                // store old values so we can check the error afterwards
                temp = v[I][j];

                v[I][j] = 2 * v[I - 1][j] - v[I - 2][j]; // vanishing gamma
				//v[I][j] = deltaS*kappa + v[I-1][j];

                // store biggest error so far;
                error = max(error, abs(temp - v[I][j]));

            }

            // for interest rate
            for (int i = 1; i < I; ++i)
            {
                // store old values so we can check the error afterwards
                temp = v[i][J];

                // set new values
                //v[k][i][J] = 0; // bond is worthless when r=infinity
                v[i][J] = v[i][J - 1]; // vanishing rho
                //v[k][i][J] = 2*v[k][i][J-1] - v[k][i][J-2];

                // store biggest error so far;
                error = max(error, abs(temp - v[i][J]));

                temp = vb[i][J];
                vb[i][J] = vb[i][J - 1];
				// rho of a call option should be zero at r = infinity

                // store biggest error so far;
                error = max(error, abs(temp - vb[i][J]));
            }

            errorOutput << error << endl;
        }
        // while post-condition: v[k] now contains approximate solution
        // which was close enough to last approximation in the iteration
        // v[k] is the result of solving the PDE implicitly using
        // the Yanenko splitting method
    }




	errorOutput.close();

    std::ofstream outputFile("convertible_prices.csv");
    // output initial convertible bond prices
	for (int i = 0; i < I+1; ++i)
	{
		for (int j =0; j < J+1; ++j)
			//outputFile << v[i][j] << ", ";
			outputFile << v[i][j] - vb[i][j] << ", ";
		outputFile << endl;
	}
	outputFile.close();
}

