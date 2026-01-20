#include <iostream>
#include <vector>
#include <cmath>

#ifndef BEAM2_SUBMODULE
#include "opengl3.h"
#else
#include "../opengl3.h"
int beam2_init(openGLframe &graphics);
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

class BeamElement3D {
private:
    double length;           // Element length (m)
    double E;                // Young's modulus (Pa)
    double G;                // Shear modulus (Pa)
    double Iyy;              // Second moment about y (m^4)
    double Izz;              // Second moment about z (m^4)
    double J;                // Polar moment / torsional constant (m^4)
    double rho;              // Density (kg/m^3)
    double area;             // Cross-sectional area (m^2)

public:
    BeamElement3D(double _length, double _E, double _G, double _Iyy, double _Izz, double _J, double _rho, double _area);

    Eigen::Matrix<double, 12, 12> getStiffnessMatrix() const;
    Eigen::Matrix<double, 12, 12> getMassMatrix() const;

};


class CantileverBeam3D {
private:
    //Model definitions:
    double length;           // Total beam length (m)
    double E;                // Young's modulus (Pa)
    double nu;               // Poisson's ratio
    double rho;              // Density (kg/m^3)
    double area;             // Cross-sectional area (m^2)
    double Iyy, Izz, J;      // section properties
    int numElements;         // Number of finite elements
    double baseOutd, IzzBase;
    double elementLength;

    Eigen::Vector3d endPointLoad; // Point load at the free end (N), 3D vector

    //Object state variables
    std::vector<BeamElement3D> elements;
    Eigen::MatrixXd globalStiffnessMatrix;
    Eigen::MatrixXd globalMassMatrix;
    Eigen::MatrixXd globalDampingMatrix;
//    Eigen::VectorXd forceVectorMag;
//    Eigen::VectorXd forceVector;
    Eigen::VectorXd u; // displacement
    Eigen::VectorXd v; // velocity

    //Simulation variables
    int totalDOFs, activeDOFs;
    Eigen::MatrixXd Kactive, Mactive, Factive, Cactive;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    Eigen::VectorXd acc;
    Eigen::LLT<Eigen::MatrixXd> lltOfLHS;
    double timeStep;
    
    // Base motion state (6 DOFs: ux, uy, uz, rotx, roty, rotz)
    Eigen::VectorXd u_base, v_base, acc_base;
    
    // Original coupling matrices (saved before boundary conditions zero them out)
    Eigen::MatrixXd K_coupling_original;  // K(active_rows, base_cols) - original before BC
    Eigen::MatrixXd M_coupling_original;   // M(active_rows, base_cols) - original before BC
    Eigen::MatrixXd K_base_to_active_original;  // K(base_rows, active_cols) - original before BC
    Eigen::MatrixXd M_base_to_active_original;   // M(base_rows, active_cols) - original before BC
    Eigen::MatrixXd K_base_base_original;  // K(base_rows, base_cols) - original before BC
    Eigen::MatrixXd M_base_base_original;   // M(base_rows, base_cols) - original before BC
    
    // Coupling matrices: forces on active DOFs due to base motion
    Eigen::MatrixXd K_coupling;  // K(active_rows, base_cols)
    Eigen::MatrixXd C_coupling;  // C(active_rows, base_cols)
    Eigen::MatrixXd M_coupling;   // M(active_rows, base_cols)
    
    // For computing reaction forces: forces on base due to active DOF motion
    Eigen::MatrixXd K_base_to_active;  // K(base_rows, active_cols)
    Eigen::MatrixXd C_base_to_active;  // C(base_rows, active_cols)
    Eigen::MatrixXd M_base_to_active;   // M(base_rows, active_cols)
    Eigen::MatrixXd K_base_base;  // K(base_rows, base_cols)
    Eigen::MatrixXd M_base_base;   // M(base_rows, base_cols)

public:
    CantileverBeam3D(double _length, double _E, double _nu, double _rho, double _area,
        int _numElements, const Eigen::Vector3d& _endPointLoad,
        double _outDia, double _inDia, double _taper);

    void solveStaticDisplacement(Eigen::VectorXd& forceVector);
    void setupTimeDomainSimulation(double _timeStep,  double dampingRatio = 0.05);
    void stepForward(double timeStep, Eigen::VectorXd& forceVector);
    void visualize(openGLframe& graphics);
    void solveFrequencyAnalysis(int numModes);
    void draw(openGLframe& graphics);
    Eigen::VectorXd applyEndpointLoad(Eigen::Vector3d endPointLoad);
    int DOF() { return totalDOFs; }
    //void setBaseOffset(Eigen::Vector3d _offset);

    void simulateTimeDomain(openGLframe& graphics, double duration, double _timeStep, double _dampingRatio = 0.05);

    void simulateTimeDomain2(openGLframe& graphics, double duration, double timeStep, double dampingRatio = 0.05);
    
    // Coupling interface for external simulation
    void setBaseState(const Eigen::VectorXd& position, 
                      const Eigen::VectorXd& velocity,
                      const Eigen::VectorXd& acceleration);
    
    // Get reaction forces at base (6 DOFs: Fx, Fy, Fz, Mx, My, Mz)
    Eigen::VectorXd getBaseReactionForces() const;
    
    // Step the beam forward with external base motion
    void stepForwardWithBaseMotion(double timeStep, 
                                   const Eigen::VectorXd& externalForces);

private:
    void showOnScreen(openGLframe& graphics, double dt= -1);
    void applyBoundaryConditions();
    void assembleGlobalMatrices();
};

CantileverBeam3D make_beam();
