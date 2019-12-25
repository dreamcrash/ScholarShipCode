#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "SeriesTest.h"



void run(const int size, const int validation){
    
    /** Matrices allocation */
    double (*TestArray) [2] = malloc(sizeof *TestArray * size);
    
    /** Dealing with fail memory allocation */
    if(!TestArray) return;
    
    double begin  = omp_get_wtime();
    
    Do(size, TestArray);
    
    const double end    = omp_get_wtime();
    
    printf("Time:%f\n",(end-begin));
    
    if(validation)
    {
        JGFvalidate(size, TestArray);
    }
    
    free(TestArray);
}

void JGFvalidate(const int size, double TestArray[2][size]){
    
    double ref[4][2] = {{2.8729524964837996, 0.0},
                       {1.1161046676147888, -1.8819691893398025},
                       {0.34429060398168704, -1.1645642623320958},
                       {0.15238898702519288, -0.8143461113044298}};

     for (int i = 0; i < 4; i++)
     {
         for (int j = 0; j < 2; j++)
         {
            double error = fabs(TestArray[j][i] - ref[i][j]);
            if (error > 1.0e-12)
            {
                printf("Validation failed for coefficient %d , %d\n",j, i);
		printf("Computed value =  %f \n",TestArray[j][i]);
		printf("Reference value = %f\n", ref[i][j]);
            }
         }
     }
  }


void Do (const int size, double TestArray [2][size]){
    
   // Calculate the fourier series. Begin by calculating A[0].
    TestArray[0][0]=TrapezoidIntegrate((double)0.0, // Lower bound.
                            (double)2.0,            // Upper bound.
                            1000,                    // # of steps.
                            (double)0.0,            // No omega*n needed.
                            0) / (double)2.0;       // 0 = term A[0].
    // Calculate the fundamental frequency.
    // ( 2 * pi ) / period...and since the period
    // is 2, omega is simply pi.

    for (int i = 1; i < size; i ++)
    {
        const double omegan = (double) 3.1415926535897932 * ((double)i);
        // Calculate A[i] terms. Note, once again, that we
        // can ignore the 2/period term outside the integral
        // since the period is 2 and the term cancels itself
        // out.

        TestArray[0][i] = TrapezoidIntegrate((double)0.0,
                          (double)2.0,
                          1000,
                          omegan,
                          1);                       // 1 = cosine term.

        // Calculate the B[i] terms.

        TestArray[1][i] = TrapezoidIntegrate((double)0.0,
                          (double)2.0,
                          1000,
                          omegan,
                          2);                       // 2 = sine term.
    }    
}


double TrapezoidIntegrate (double x0,     // Lower bound.
                        double x1,                // Upper bound.
                        int nsteps,               // # of steps.
                        double omegan,            // omega * n.
                        int select)               // Term type.
{
    double x;               // Independent variable.
    double dx;              // Step size.
    double rvalue;          // Return value.

    // Initialize independent variable.

    x = x0;

    // Calculate stepsize.

    dx = (x1 - x0) / (double)nsteps;

    // Initialize the return value.

    rvalue = thefunction(x0, omegan, select) / (double)2.0;

    // Compute the other terms of the integral.

    if (nsteps != 1)
    {
            --nsteps;               // Already done 1 step.
            while (--nsteps > 0)
            {
                    x += dx;
                    rvalue += thefunction(x, omegan, select);
            }
    }

    // Finish computation.

    rvalue=(rvalue + thefunction(x1,omegan,select) / (double)2.0) * dx;
    return(rvalue);
}

/*
* thefunction
*
* This routine selects the function to be used in the Trapezoid
* integration. x is the independent variable, omegan is omega * n,
* and select chooses which of the sine/cosine functions
* are used. Note the special case for select=0.
*/

double thefunction(double x,      // Independent variable.
                double omegan,              // Omega * term.
                int select)                 // Choose type.
{

    // Use select to pick which function we call.

    switch(select)
    {
        case 0: return(pow(x+(double)1.0,x));

        case 1: return(pow(x+(double)1.0,x) * cos(omegan*x));

        case 2: return(pow(x+(double)1.0,x) * sin(omegan*x));
    }

    // We should never reach this point, but the following
    // keeps compilers from issuing a warning message.

    return (0.0);
}