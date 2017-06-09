//
//  updPosQ.c
//  FreshStart
//
//  Created by Gideon Haber on 6/8/17.
//  Copyright © 2017 Summer17. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "langHeader2.h"

//UPDATE THE POSITION OF A PARTICLE
double updPosQ(double pos, double vel, double C) {
    double newPos = pos + TIMESTEP * vel + C;
    return newPos;
}
