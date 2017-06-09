//
//  langHeader2.h
//  FreshStart
//
//  Created by Gideon Haber on 6/8/17.
//  Copyright © 2017 Summer17. All rights reserved.
//

#ifndef langHeader2_h
#define langHeader2_h
#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

const double BOLTZMANN = 0.05;
const double TEMP = 200;
const int NUM_TRAJ = 100;
const int NUM_PART = 1;
const int HBAR = 1;
const double SPRING_F = NUM_PART * BOLTZMANN * TEMP / HBAR;
const double OMEGA = SPRING_F / NUM_PART;
const double MASS = 1;
const double GAMMA = 0.1;
const double TIMESTEP = .00001;
const double SIMULENGTH = 2;

// GENERATE UNIFORMLY DISTRIBUTED RANDOM
double genRandom(){
    double rNum = (double)rand() / (double)((unsigned)RAND_MAX + 1);
    return rNum;
}


// RETURN POSITION OF THE NEXT PARTICLE IN RING
// tested and functional :)
double nextPos(double* column, int particle) {
    double npos;
    if (particle == (NUM_PART - 1)) {
        npos = *column;
    }
    else {
        npos = *(column + particle + 1);
    }
    return npos;
}

// RETURN POSITION OF THE PREVIOUS PARTICLE IN RING
double prevPos(double* column, int particle) {
    double ppos;
    if (particle == 0) {
        ppos = *(column + NUM_PART - 1);
    }
    else {
        ppos = *(column + particle - 1);
    }
    return ppos;
}

// GENERATE GAUSSIAN DISTRIBUTED VARIABLES USING BOX-MUELLER FORMULA
double * genGauss() {
    static double numbers[2];
    double U = genRandom();
    double V = genRandom();
    double X = sqrt(-2 * log(U)) * cos(2 * M_PI * V);
    double Y = sqrt(-2 * log(U)) * sin(2 * M_PI * V);
    numbers[0] = X, numbers[1] = Y;
    return numbers;
}


double calcOmegaL(int l){
    double om = 2.0 * SPRING_F * sin(l * M_PI / NUM_PART);
    return om;
}

// GENERATE INITIAL MOMENTA AND POSITIONS (NOT SURE IF CORRECT...)
double* genIntialCond(int l) {
    // set up space for array whose pointer will be returned
    static double numbers[2];
    
    // generate random numbers
    double * rNump = genGauss();
    double X = *rNump;
    double Y = *(rNump + 1);
    
    // find value of omega to initalize
    double om = calcOmegaL(l);
    
    // convert mean/std Dev of random pos/vel
    Y = Y * sqrt(BOLTZMANN * TEMP * NUM_PART);
    if (om == 0)
        X = 0;
    else
        X = X * sqrt(BOLTZMANN * TEMP / pow(om, 2));
    numbers[0] = X, numbers[1] = Y;
    
    // return pointer to random intial pos/vel
    return numbers;
}

// RETURN POINTER TO THE COLUMN OF AN ARRAY
double* retColumn(double matrix[NUM_PART][NUM_TRAJ], int col){
    static double column[NUM_PART];
    for (int j = 0; j < NUM_PART; j++){
        column[j] = matrix[j][col];
    }
    return column;
}

// RETURN POINTER TO THE ROW OF AN ARRAY
double* retRow(double matrix[NUM_PART][NUM_TRAJ], int row){
    static double rowHold[NUM_TRAJ];
    for (int j = 0; j < NUM_TRAJ; j++){
        rowHold[j] = matrix[row][j];
    }
    return rowHold;
}

// FUNCTION FOR CALCULATION FORCE
double calcForce(double* column, int j) {
    double pos = *(column + j);
    double npos = nextPos(column, j);
    double ppos = prevPos(column, j);
    double force = 1.0 / MASS * (-pow(SPRING_F, 2) *
                    (2.0 * pos - npos - ppos) - pow(OMEGA, 2) * pos);
    force = 0;
    return force;
}

//FUNCTION FOR CALCULATING C VARIABLE
double * calcC(double force, double vel){
    // create constant
    static double cPlus[2];
    double sigma = sqrt(2.0 * BOLTZMANN * TEMP * GAMMA / MASS);
    
    // generate two random variables
    double* rVar = genGauss();
    double xi = *rVar;
    double theta = *(rVar + 1);
    
    // begin calculation of C (break up since otherwise too long)
    double C1 = 1.0 / 2 * pow(TIMESTEP, 2) * (force - GAMMA * vel);
    double C2 = sigma * pow(TIMESTEP, 1.5) * (0.5 * xi + 1.0 / (2.0 * sqrt(3)) * theta);
    cPlus[0] = C1 + C2, cPlus[1] = xi;
    return cPlus;
}

//FUNCTION FOR CALCULATING THE AVERAGE ENERGY OF AN ARRAY
double calcAveE(double* velocity, double* position) {
    double holder = 0;
    for (int i = 0; i < NUM_TRAJ; i++) {
        double u = 1.0 / 2 * pow(*(velocity + i), 2);
//        holder = holder + u + 1.0 / 2 * pow(OMEGA, 2) * pow(*(position + i), 2);
        holder = holder + u;
    }
    double avg = holder / NUM_TRAJ;
    return avg;
}

double updPosQ(double pos, double vel, double C);
double updVelQ(double cForce, double pForce, double vel, double cPlus[]);

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


//UPDATE THE POSITION OF A PARTICLE
double updPosQ(double pos, double vel, double C) {
    double newPos = pos + TIMESTEP * vel + C;
    return newPos;
}


#endif /* langHeader2_h */
