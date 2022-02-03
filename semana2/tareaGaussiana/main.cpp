//
//  main.cpp
//  tutorialcpp
//
//  Created by Carlos Alberto Salas Salas on 7/7/21.
//  Copyright © 2021 Carlos Alberto Salas Salas. All rights reserved.
//

#include <iostream>
#include <cmath>

using namespace std;


double gaussian(double x){
    return exp(-x*x);
}

double centralDerivative(double x, double h = 1E-6){
    return (gaussian(x + h) - gaussian(x - h))/(2*h);
}

double obtainDerivative(double xMin, double xMax, double step = 0.01){
    
    int stepsInt = round((xMax - xMin)/step) + 1;
    
    double xLinspace [stepsInt];
    double yDerivativeValues [stepsInt];
    
    for (int x = 0; x < stepsInt; x++) {
        xLinspace[x] = xMin + x*step;
        yDerivativeValues[x] = centralDerivative(xLinspace[x]);
        cout << xLinspace[x] << "," <<  yDerivativeValues[x] << '\n';
    }
    
    
    
    return 0.;
    
}



int main(int argc, const char * argv[]) {
    
    double xMin = -20.;
    double xMax = 20.;
    
    cout << "x,yPrime" << '\n';
    
    obtainDerivative(xMin, xMax);
    

    

    return 0;

}

