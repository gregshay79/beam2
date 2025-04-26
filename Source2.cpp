//#include <windows.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "opengl3.h"
#include "heatmap.h"


#define M_PI 3.141592653589793238462643383 

mygraphics graphics;

// Function to calculate moment of inertia for a hollow circular tube
double calculateHollowTubeMomentOfInertia(double outerDiameter, double innerDiameter) {
    double outerRadius = outerDiameter / 2.0;
    double innerRadius = innerDiameter / 2.0;
    return M_PI * (pow(outerRadius, 4) - pow(innerRadius, 4)) / 4.0;
}

// Function to calculate cross-sectional area for a hollow circular tube
double calculateHollowTubeArea(double outerDiameter, double innerDiameter) {
    double outerRadius = outerDiameter / 2.0;
    double innerRadius = innerDiameter / 2.0;
    return M_PI * (pow(outerRadius, 2) - pow(innerRadius, 2));
}

class BeamElement {
private:
    double length;           // Element length (m)
    double E;                // Young's modulus (Pa)
    double I;                // Second moment of area (m^4)
    double rho;              // Density (kg/m^3)
    double area;             // Cross-sectional area (m^2)

public:
    BeamElement(double _length, double _E, double _I, double _rho, double _area)
        : length(_length), E(_E), I(_I), rho(_rho), area(_area) {
    }

    // Local stiffness matrix for beam element (4x4)
    Eigen::Matrix4d getStiffnessMatrix() const {
        Eigen::Matrix4d K = Eigen::Matrix4d::Zero();

        double factor = E * I / (length * length * length);

        K(0, 0) = 12.0;
        K(0, 1) = 6.0 * length;
        K(0, 2) = -12.0;
        K(0, 3) = 6.0 * length;

        K(1, 0) = 6.0 * length;
        K(1, 1) = 4.0 * length * length;
        K(1, 2) = -6.0 * length;
        K(1, 3) = 2.0 * length * length;

        K(2, 0) = -12.0;
        K(2, 1) = -6.0 * length;
        K(2, 2) = 12.0;
        K(2, 3) = -6.0 * length;

        K(3, 0) = 6.0 * length;
        K(3, 1) = 2.0 * length * length;
        K(3, 2) = -6.0 * length;
        K(3, 3) = 4.0 * length * length;

        return factor * K;
    }

    // Local mass matrix for beam element (4x4) - consistent mass matrix
    Eigen::Matrix4d getMassMatrix() const {
        Eigen::Matrix4d M = Eigen::Matrix4d::Zero();

        double factor = rho * area * length / 420.0;

        M(0, 0) = 156.0;
        M(0, 1) = 22.0 * length;
        M(0, 2) = 54.0;
        M(0, 3) = -13.0 * length;

        M(1, 0) = 22.0 * length;
        M(1, 1) = 4.0 * length * length;
        M(1, 2) = 13.0 * length;
        M(1, 3) = -3.0 * length * length;

        M(2, 0) = 54.0;
        M(2, 1) = 13.0 * length;
        M(2, 2) = 156.0;
        M(2, 3) = -22.0 * length;

        M(3, 0) = -13.0 * length;
        M(3, 1) = -3.0 * length * length;
        M(3, 2) = -22.0 * length;
        M(3, 3) = 4.0 * length * length;

        return factor * M;
    }
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
    CantileverBeam(double _length, double _E, double _rho, double _area, int _numElements, double _pointLoad, double _outDia, double _inDia, double _taper)
        : length(_length), E(_E),/* I(_I),*/ rho(_rho), area(_area), numElements(_numElements), pointLoad(_pointLoad) {

        // Initialize elements
        double elementLength = length / numElements;
        for (int i = 0; i < numElements; ++i) {

            double taper_factor = 1.0 - _taper * i / numElements;  // taper facter from 1.0 down to (1 - _taper).
            double I = calculateHollowTubeMomentOfInertia(_outDia * taper_factor, _inDia * taper_factor);
            double A = calculateHollowTubeArea(_outDia * taper_factor, _inDia * taper_factor);
            elements.emplace_back(elementLength, E, I, rho, A);
        }

        // Size of the global matrices: 2 DOFs per node (deflection and rotation)
        int numDOFs = 2 * (numElements + 1);
        globalStiffnessMatrix = Eigen::MatrixXd::Zero(numDOFs, numDOFs);
        globalMassMatrix = Eigen::MatrixXd::Zero(numDOFs, numDOFs);
        forceVector = Eigen::VectorXd::Zero(numDOFs);

        assembleGlobalMatrices();
        applyBoundaryConditions();
        applyLoads();
    }

    void assembleGlobalMatrices() {
        for (int i = 0; i < numElements; ++i) {
            Eigen::Matrix4d localK = elements[i].getStiffnessMatrix();
            Eigen::Matrix4d localM = elements[i].getMassMatrix();

            // Map local DOFs to global DOFs
            int startDOF = 2 * i;  // First DOF of the element in global numbering

            // Add local matrices to global matrices
            for (int r = 0; r < 4; ++r) {
                for (int c = 0; c < 4; ++c) {
                    int globalR = startDOF + (r / 2) * 2 + (r % 2);
                    int globalC = startDOF + (c / 2) * 2 + (c % 2);

                    globalStiffnessMatrix(globalR, globalC) += localK(r, c);
                    globalMassMatrix(globalR, globalC) += localM(r, c);
                }
            }
        }
    }

    void applyBoundaryConditions() {
        // For a cantilever beam, fix the DOFs at the first node (x=0)
        // We use a penalty method by setting very large values in the stiffness matrix
        double penaltyFactor = 1.0e15;
        globalStiffnessMatrix(0, 0) = penaltyFactor;  // Fix deflection
        globalStiffnessMatrix(1, 1) = penaltyFactor;  // Fix rotation

        // Zero out coupling terms
        for (int i = 0; i < globalStiffnessMatrix.rows(); ++i) {
            if (i != 0) globalStiffnessMatrix(0, i) = 0.0;
            if (i != 1) globalStiffnessMatrix(1, i) = 0.0;
            if (i != 0) globalStiffnessMatrix(i, 0) = 0.0;
            if (i != 1) globalStiffnessMatrix(i, 1) = 0.0;
        }
    }

    void applyLoads() {
        // Apply point load at the free end (last node, vertical DOF)
        int lastNodeDOF = 2 * numElements;
        forceVector(lastNodeDOF) = pointLoad;
    }

    void solveStaticDisplacement() {
        // Solve for static displacements: K * u = F
        // First, extract the active DOFs (excluding fixed DOFs at first node)
        int numDOFs = 2 * (numElements + 1) - 2;
        Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(numDOFs, numDOFs);
        Eigen::VectorXd Factive = forceVector.tail(numDOFs);

        // Solve the system
        Eigen::VectorXd displacements = Kactive.colPivHouseholderQr().solve(Factive);

        // Output results
        std::cout << "Static Displacement Analysis Results:" << std::endl;
        std::cout << "-----------------------------------" << std::endl;

        // Full displacement vector including fixed DOFs
        Eigen::VectorXd fullDisplacements = Eigen::VectorXd::Zero(2 * (numElements + 1));
        fullDisplacements.tail(numDOFs) = displacements;

        for (int i = 0; i <= numElements; ++i) {
            double position = i * (length / numElements);
            std::cout << "Node " << i << " (x = " << position << " m):" << std::endl;
            std::cout << "  Deflection: " << fullDisplacements(2 * i) << " m" << std::endl;
            std::cout << "  Rotation: " << fullDisplacements(2 * i + 1) << " rad" << std::endl;
            std::cout << std::endl;
        }
    }

    void solveFrequencyAnalysis(int numModes) {
        // Extract active DOFs (excluding fixed DOFs)
        int numDOFs = 2 * (numElements + 1) - 2;
        Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(numDOFs, numDOFs);
        Eigen::MatrixXd Mactive = globalMassMatrix.bottomRightCorner(numDOFs, numDOFs);

        // Solve the generalized eigenvalue problem: K * v = λ * M * v
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(Kactive, Mactive);

        // Output natural frequencies
        std::cout << "Natural Frequencies:" << std::endl;
        std::cout << "-------------------" << std::endl;

        for (int i = 0; i < std::min(numModes, numDOFs); ++i) {
            double naturalFreq = std::sqrt(solver.eigenvalues()(i)) / (2.0 * M_PI); // in Hz
            std::cout << "Mode " << (i + 1) << ": " << naturalFreq << " Hz" << std::endl;
        }

        // Output mode shapes if desired
        std::cout << "\nMode Shapes (normalized):" << std::endl;
        std::cout << "----------------------" << std::endl;

        for (int mode = 0; mode < std::min(numModes, numDOFs); ++mode) {
            Eigen::VectorXd modeShape = solver.eigenvectors().col(mode);

            // Normalize the mode shape
            modeShape /= modeShape.norm();

            std::cout << "Mode " << (mode + 1) << ":" << std::endl;

            // Create full mode shape including fixed DOFs
            Eigen::VectorXd fullModeShape = Eigen::VectorXd::Zero(2 * (numElements + 1));
            fullModeShape.tail(numDOFs) = modeShape;

            for (int i = 0; i <= numElements; ++i) {
                double position = i * (length / numElements);
                std::cout << "  Node " << i << " (x = " << position << " m):" << std::endl;
                std::cout << "    Deflection: " << fullModeShape(2 * i) << std::endl;
                std::cout << "    Rotation: " << fullModeShape(2 * i + 1) << std::endl;
            }
            std::cout << std::endl;
        }
    }

    // Compute the derivative for the RK4 integration
    void computeStateDerivative(
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& v,
        const Eigen::MatrixXd& invM,
        const Eigen::MatrixXd& K,
        const Eigen::MatrixXd& C,
        const Eigen::VectorXd& F,
        Eigen::VectorXd& dudt,
        Eigen::VectorXd& dvdt) {

        // dudt = v
        dudt = v;

        // dvdt = M^-1 * (F - K*u - C*v)
        dvdt = invM * (F - K * u - C * v);
    }

    //void simulateTimeDomainRK4(double duration, double timeStep, double dampingRatio = 0.05) {
    //    // Extract active DOFs (excluding fixed DOFs)
    //    int numDOFs = 2 * (numElements + 1) - 2;
    //    Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(numDOFs, numDOFs);
    //    Eigen::MatrixXd Mactive = globalMassMatrix.bottomRightCorner(numDOFs, numDOFs);
    //    Eigen::VectorXd Factive = forceVector.tail(numDOFs);

    //    // Compute Rayleigh damping matrix: C = alpha*M + beta*K
    //    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(Kactive, Mactive);
    //    double omega1 = std::sqrt(solver.eigenvalues()(0));
    //    double omega2 = std::sqrt(solver.eigenvalues()(1));

    //    // Solve for alpha and beta to get desired damping at two frequencies
    //    double alpha = dampingRatio * 2.0 * omega1 * omega2 / (omega1 + omega2);
    //    double beta = dampingRatio * 2.0 / (omega1 + omega2);

    //    Eigen::MatrixXd Cactive = alpha * Mactive + beta * Kactive;

    //    // Pre-compute the inverse of mass matrix for efficiency
    //    Eigen::MatrixXd invMactive = Mactive.inverse();

    //    // Initial conditions (zero displacement and velocity)
    //    Eigen::VectorXd u = Eigen::VectorXd::Zero(numDOFs);
    //    Eigen::VectorXd v = Eigen::VectorXd::Zero(numDOFs);

    //    std::cout << "Time-Domain Simulation (RK4):" << std::endl;
    //    std::cout << "-------------------------" << std::endl;
    //    std::cout << "Time (s), Tip Deflection (m)" << std::endl;

    //    // Time stepping with 4th order Runge-Kutta
    //    int numTimeSteps = static_cast<int>(duration / timeStep);

    //    // Temporary variables for the RK4 method
    //    Eigen::VectorXd k1u, k1v, k2u, k2v, k3u, k3v, k4u, k4v;
    //    k1u.resize(numDOFs); k1v.resize(numDOFs);
    //    k2u.resize(numDOFs); k2v.resize(numDOFs);
    //    k3u.resize(numDOFs); k3v.resize(numDOFs);
    //    k4u.resize(numDOFs); k4v.resize(numDOFs);

    //    for (int step = 0; step <= numTimeSteps; ++step) {
    //        double time = step * timeStep;

    //        // Print tip deflection (last displacement DOF)
    //        std::cout << time << ", " << u(numDOFs - 2) << std::endl;

    //        // RK4 integration step



    //        // Step 1: Evaluate derivatives at current time
    //        computeStateDerivative(u, v, invMactive, Kactive, Cactive, Factive, k1u, k1v);

    //        // Step 2: Evaluate derivatives at time + h/2 using step 1 derivatives
    //        computeStateDerivative(
    //            u + 0.5 * timeStep * k1u,
    //            v + 0.5 * timeStep * k1v,
    //            invMactive, Kactive, Cactive, Factive, k2u, k2v);

    //        // Step 3: Evaluate derivatives at time + h/2 using step 2 derivatives
    //        computeStateDerivative(
    //            u + 0.5 * timeStep * k2u,
    //            v + 0.5 * timeStep * k2v,
    //            invMactive, Kactive, Cactive, Factive, k3u, k3v);

    //        // Step 4: Evaluate derivatives at time + h using step 3 derivatives
    //        computeStateDerivative(
    //            u + timeStep * k3u,
    //            v + timeStep * k3v,
    //            invMactive, Kactive, Cactive, Factive, k4u, k4v);

    //        // Update displacements and velocities using weighted average of the four derivatives
    //        u += (timeStep / 6.0) * (k1u + 2.0 * k2u + 2.0 * k3u + k4u);
    //        v += (timeStep / 6.0) * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);

    //        // Optional: Apply time-varying force here if needed
    //    }
    //}

    void simulateTimeDomain(double duration, double timeStep, double dampingRatio = 0.05) {

        // Vector to store line segments
        std::vector<LineSegment> lineSegments;

        // Extract active DOFs (excluding fixed DOFs)
        int numDOFs = 2 * (numElements + 1) - 2;
        Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(numDOFs, numDOFs);
        Eigen::MatrixXd Mactive = globalMassMatrix.bottomRightCorner(numDOFs, numDOFs);
        Eigen::VectorXd Factive = forceVector.tail(numDOFs);

        // Compute Rayleigh damping matrix: C = alpha*M + beta*K
        // We'll use modal damping for simplicity with the same ratio for all modes
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(Kactive, Mactive);
        double omega1 = std::sqrt(solver.eigenvalues()(0));
        double omega2 = std::sqrt(solver.eigenvalues()(1));

        // Solve for alpha and beta to get desired damping at two frequencies
        double alpha = dampingRatio * 2.0 * omega1 * omega2 / (omega1 + omega2);
        double beta = dampingRatio * 2.0 / (omega1 + omega2);

        Eigen::MatrixXd Cactive = alpha * Mactive + beta * Kactive;

        // Initial conditions (zero displacement and velocity)
        Eigen::VectorXd u = Eigen::VectorXd::Zero(numDOFs);
        Eigen::VectorXd v = Eigen::VectorXd::Zero(numDOFs);
        Eigen::VectorXd a = Mactive.colPivHouseholderQr().solve(Factive - Kactive * u - Cactive * v);

        // Newmark-beta time integration parameters
        double gamma = 0.5;
        double beta_nb = 0.25;  // Not to be confused with Rayleigh damping beta

        // Constants for Newmark-beta method
        Eigen::MatrixXd LHS = Mactive + gamma * timeStep * Cactive + beta_nb * timeStep * timeStep * Kactive;
        Eigen::LLT<Eigen::MatrixXd> lltOfLHS(LHS);  // For efficient repeated solves

        std::cout << "Time-Domain Simulation:" << std::endl;
        std::cout << "---------------------" << std::endl;
        std::cout << "Time (s), Tip Deflection (m)" << std::endl;
        float maxbend = 0;
        int numTimeSteps = static_cast<int>(duration / timeStep);
        for (int step = 0; step <= numTimeSteps; ++step) {
            double time = step * timeStep;

            // Print tip deflection (last displacement DOF)
            std::cout << time << ", " << u(numDOFs - 2) << std::endl;

            lineSegments.clear();
            float sc = 1.9/numDOFs;
            for (int i = 0; i < numDOFs - 4; i+=2) {
                float x1, y1, x2, y2;
                x1 = sc * i - 0.9;
                y1 = -1.0 * u(i);
                x2 = sc * (i + 2) - 0.9;
                y2 = -1.0 * u(i+2);
                float bend = 2000.0/1.6 * fabs(u(i + 1)-u(i+3));
                if (bend > maxbend) maxbend = bend;
                RGB color = valueToHeatmapColor(bend);
                lineSegments.emplace_back(x1,y1,x2,y2, color.r/255.0, color.g/255.0, color.b/255.0); 
            }

            graphics.draw(lineSegments);
            //Sleep(100);

            // Predictor step for Newmark-beta method
            Eigen::VectorXd uPred = u + timeStep * v + timeStep * timeStep * (0.5 - beta_nb) * a;
            Eigen::VectorXd vPred = v + timeStep * (1.0 - gamma) * a;

            // RHS for the corrector step
            Eigen::VectorXd RHS = Factive - Kactive * uPred - Cactive * vPred;

            // Corrector step
            Eigen::VectorXd deltaA = lltOfLHS.solve(RHS);

            // Update acceleration, velocity, and displacement
            a = deltaA;
            v = vPred + gamma * timeStep * a;
            u = uPred + beta_nb * timeStep * timeStep * a;

            // Optional: Apply time-varying force here if needed
        }
        std::cout << "max bend = " << maxbend;
    }


};



int main() {
    
    graphics.setupGL();

    // Hollow aluminum tube parameters
    double length = 12.8;                 // Length (m)
    double thickness = .125 / 39.37;
    double outerDiameter = 8/39.37;         // Outer diameter (m)
    double innerDiameter = outerDiameter - 2*thickness;        // Inner diameter (m)
    double E = 69e9;                     // Young's modulus for aluminum (Pa)
    double rho = 2700.0;                 // Density of aluminum (kg/m^3)
    double taper = 0.5;

    // Calculate properties for a hollow circular tube
    double area = calculateHollowTubeArea(outerDiameter, innerDiameter);
    double I = calculateHollowTubeMomentOfInertia(outerDiameter, innerDiameter);

    int numElements = 40;               // Number of finite elements
    double pointLoad = 10.0*4.448;           // Point load at the free end (N)

    std::cout << "Finite Element Analysis of a Hollow Aluminum Tube Cantilever" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "Beam length: " << length << " m" << std::endl;
    std::cout << "Outer diameter: " << outerDiameter << " m" << std::endl;
    std::cout << "Inner diameter: " << innerDiameter << " m" << std::endl;
    std::cout << "Wall thickness: " << (outerDiameter - innerDiameter) / 2.0 << " m" << std::endl;
    std::cout << "Cross-sectional area: " << area << " m^2" << std::endl;
    std::cout << "Moment of inertia: " << I << " m^4" << std::endl;
    std::cout << "Young's modulus: " << E << " Pa" << std::endl;
    std::cout << "Density: " << rho << " kg/m^3" << std::endl;
    std::cout << "Point load at free end: " << pointLoad << " N" << std::endl;
    std::cout << "Number of elements: " << numElements << std::endl << std::endl;

    // Create and analyze the beam
    CantileverBeam beam(length, E, rho, area, numElements, pointLoad, outerDiameter, innerDiameter, taper);

    // Static analysis
    beam.solveStaticDisplacement();

    std::cout << std::endl;

    // Modal analysis (first 3 modes)
    //beam.solveFrequencyAnalysis(3);

    std::cout << std::endl;

    // Time domain simulation for 0.5 seconds with 0.01s time step using RK4
    //beam.simulateTimeDomainRK4(0.5, 0.01);

    // Time domain simulation for 2 seconds with 0.01s time step
    beam.simulateTimeDomain(10.0, 0.05);

    //for (int i = 0; i < 1000; i++) {
    //    //Sleep(100);
    //    for (auto& line : lineSegments) {
    //        line.x1 *= .99;
    //        line.x2 *= .99;
    //        line.y1 *= .99;
    //        line.y2 *= .99;
    //    }
    //    graphics.draw(lineSegments);

    //    glfwPollEvents();
    //    if (glfwWindowShouldClose(graphics.window))
    //        break;
    //}

    graphics.waitForCompletion();

    graphics.closeGL();

    return 0;
}