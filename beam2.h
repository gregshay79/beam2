#include <iostream>
#include <vector>
#include <cmath>

#ifndef BEAM2_SUBMODULE
#include "opengl3.h"
#else
#include "../opengl3.h"
#endif

double calculateHollowTubeMomentOfInertia(double outerDiameter, double innerDiameter);
double calculateHollowTubeArea(double outerDiameter, double innerDiameter);

class BeamElement {
private:
    double length;           // Element length (m)
    double E;                // Young's modulus (Pa)
    double I;                // Second moment of area (m^4)
    double rho;              // Density (kg/m^3)
    double area;             // Cross-sectional area (m^2)

public:
    BeamElement(double _length, double _E, double _I, double _rho, double _area);
    Eigen::Matrix4d getStiffnessMatrix() const;
    Eigen::Matrix4d getMassMatrix() const;
};

class CantileverBeam {
private:
    double length;           // Total beam length (m)
    double E;                // Young's modulus (Pa)
    //double I;                // Second moment of area (m^4)
    double rho;              // Density (kg/m^3)
    double area;             // Cross-sectional area (m^2)
    int numElements;         // Number of finite elements
    double pointLoad;        // Point load at the free end (N)

    std::vector<BeamElement> elements;
    Eigen::MatrixXd globalStiffnessMatrix;
    Eigen::MatrixXd globalMassMatrix;
    Eigen::MatrixXd globalDampingMatrix;
    Eigen::VectorXd forceVector;

    // Pre-factorized mass matrix inverse for more efficient time integration
    Eigen::LLT<Eigen::MatrixXd> massMatrixSolver;

public:
    CantileverBeam(double _length, double _E, double _rho, double _area, int _numElements, double _pointLoad, double _outDia, double _inDia, double _taper);
    
    void simulateTimeDomain(openGLframe &_graphics, double duration, double timeStep, double dampingRatio = 0.05);


public:
    void assembleGlobalMatrices();
    void applyBoundaryConditions();
    void applyLoads();
    void solveStaticDisplacement();
    void solveFrequencyAnalysis(int numModes);


};
