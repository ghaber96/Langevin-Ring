//
//  main.c
//  FreshStart
//
//  Created by Gideon Haber on 6/8/17.
//  Copyright © 2017 Summer17. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "langHeader2.h"

int main() {
    srand((unsigned)time(NULL));\
    
    // open file to print out data
    FILE *out_file = fopen("/Users/gideon/Desktop/debugPlease", "w");
    if (out_file == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1); // must include stdlib.h
    }
    
    
    // set up arrays to hold particles and trajectories
    double posArray[NUM_PART][NUM_TRAJ] = {0};
    double velArray[NUM_PART][NUM_TRAJ] = {0};
    double cForce[NUM_PART][NUM_TRAJ] = {0};

    // initialize arrays
    for (int h = 0; h < NUM_TRAJ; h++) {
        for (int l = 0; l < NUM_PART; l++) {
            double* cond = genIntialCond(l);
            posArray[l][h] = *cond;
            velArray[l][h] = *(cond + 1);
        }
        double* columnI = retColumn(posArray, h);
        for (int j = 0; j < NUM_PART; j++) {
            cForce[j][h] = calcForce(columnI, j);
        }
    }
    // go through times
    for(double time = 0; time < SIMULENGTH; time += TIMESTEP) {
        // go through trajectors
        for (int trajN = 0; trajN < NUM_TRAJ; trajN++) {
            // get a snap of all positions of a given trajectory to calc force
            double* column = retColumn(posArray, trajN);
            double cStore[NUM_PART][2];
            //update the positions of all particles of a given trajectory
            for (int j = 0; j < NUM_PART; j++) {
                double velHolder = velArray[j][trajN];
                double* C = calcC(cForce[j][trajN], velHolder);
                cStore[j][0] = *C, cStore[j][1] = *(C + 1);
                posArray[j][trajN] = updPosQ(posArray[j][trajN], velHolder, *C);
            }
            // get snap of all positions at new trajectory to get f(r + dt)
            column = retColumn(posArray, trajN);
            
            // update velocity
            for (int j = 0; j < NUM_PART; j++) {
                double velHolder = velArray[j][trajN];
                double nForce = calcForce(column, j);
                double cOfInt[2] = {cStore[j][0], cStore[j][1]};
                velArray[j][trajN] = updVelQ(nForce, cForce[j][trajN], velHolder, cOfInt);
                cForce[j][trajN] = nForce;
                
            }

            
        }
    }
    
    for (int i = 0; i < NUM_PART; i++) {
        double* handyVel = retRow(velArray, i);
        double* handyPos = retRow(posArray, i);
        double avgE = calcAveE(handyVel, handyPos);
        printf("%f\n", avgE);
    }

    

//    double test[NUM_TRAJ] = {1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//    double* test1 = test;
//    printf("%f", calcAveE(test1, test1));
//    
    return 0;
}
