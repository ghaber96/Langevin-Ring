//
//  updVelQ.c
//  FreshStart
//
//  Created by Gideon Haber on 6/8/17.
//  Copyright © 2017 Summer17. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "langHeader2.h"

// UPDATE THE VELOCITY OF A PARTICLE
double updVelQ(double cForce, double pForce, double vel, double cPlus[]){
    
    // define constants and get random varaibles
    double sigma = pow(2.0 * BOLTZMANN * TEMP * GAMMA / MASS, 0.5);
    
    // calculate new velocity (newVel)
    // cplus[1] = xi
    double newVel1 = vel + TIMESTEP / 2.0 * (cForce + pForce) - TIMESTEP * GAMMA * vel;
    double newVel2 = sigma * sqrt(TIMESTEP) * cPlus[1] - GAMMA * cPlus[0];
    double newVel = newVel1 + newVel2;
    return newVel;
}
