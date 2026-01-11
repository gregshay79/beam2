#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <cstdint>
#include <thread>
#include <chrono>
#include "beam2.h"

#ifndef BEAM2_SUBMODULE
#include "opengl3.h"
#include "heatmap.h"
#else
#include "../opengl3.h"
#include "../heatmap.h"
#endif

#define M_PI 3.141592653589793238462643383


double thickness = 0.525 / 39.37;              // 1/4" -> m


// Function to calculate second moment of area for a hollow circular tube
double calculateHollowTubeMomentOfInertia(double outerDiameter, double innerDiameter) {
    return M_PI * (pow(outerDiameter, 4) - pow(innerDiameter, 4)) / 64.0;
}

// Function to calculate cross-sectional area for a hollow circular tube
double calculateHollowTubeArea(double outerDiameter, double innerDiameter) {
    double outerRadius = outerDiameter / 2.0;
    double innerRadius = innerDiameter / 2.0;
    return M_PI * (pow(outerRadius, 2) - pow(innerRadius, 2));
}

//class BeamElement3D {
//private:
//    double length;           // Element length (m)
//    double E;                // Young's modulus (Pa)
//    double G;                // Shear modulus (Pa)
//    double Iyy;              // Second moment about y (m^4)
//    double Izz;              // Second moment about z (m^4)
//    double J;                // Polar moment / torsional constant (m^4)
//    double rho;              // Density (kg/m^3)
//    double area;             // Cross-sectional area (m^2)
//
//public:
    BeamElement3D::BeamElement3D(double _length, double _E, double _G, double _Iyy, double _Izz, double _J, double _rho, double _area)
        : length(_length), E(_E), G(_G), Iyy(_Iyy), Izz(_Izz), J(_J), rho(_rho), area(_area) {
    }

    // Local stiffness matrix for 3D Euler-Bernoulli beam element (12x12)
    Eigen::Matrix<double, 12, 12> BeamElement3D::getStiffnessMatrix() const {
        Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
        double L = length;
        double L2 = L * L;
        double L3 = L2 * L;

        // Axial
        double EA_L = E * area / L;
        K(0, 0) = EA_L; K(0, 6) = -EA_L;
        K(6, 0) = -EA_L; K(6, 6) = EA_L;

        // Torsion
        double GJ_L = G * J / L;
        K(3, 3) = GJ_L; K(3, 9) = -GJ_L;
        K(9, 3) = -GJ_L; K(9, 9) = GJ_L;

        // Bending about z  -> transverse displacement in y (DOFs 1,5,7,11)
        double kyz = E * Izz / L3;
        int idx_y[4] = { 1, 5, 7, 11 };
        double Kyz_local[4][4] = {
            {12.0, 6.0 * L, -12.0, 6.0 * L},
            {6.0 * L, 4.0 * L2, -6.0 * L, 2.0 * L2},
            {-12.0, -6.0 * L, 12.0, -6.0 * L},
            {6.0 * L, 2.0 * L2, -6.0 * L, 4.0 * L2}
        };
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                K(idx_y[r], idx_y[c]) += kyz * Kyz_local[r][c];

        // Bending about y  -> transverse displacement in z (DOFs 2,4,8,10)
        double kzy = E * Iyy / L3;
        int idx_z[4] = { 2, 4, 8, 10 };
        double Kzy_local[4][4] = {
            {12.0, -6.0 * L, -12.0, -6.0 * L},
            {-6.0 * L, 4.0 * L2, 6.0 * L, 2.0 * L2},
            {-12.0, 6.0 * L, 12.0, 6.0 * L},
            {-6.0 * L, 2.0 * L2, 6.0 * L, 4.0 * L2}
        };
        // Note: signs arranged to match standard element stiffness orientation
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                K(idx_z[r], idx_z[c]) += kzy * Kzy_local[r][c];

        return K;
    }

    // Simple lumped mass matrix for beam element (12x12)
    // This is a pragmatic choice to keep the implementation compact; replace with a consistent mass matrix if higher accuracy is required.
    Eigen::Matrix<double, 12, 12> BeamElement3D::getMassMatrix() const {
        Eigen::Matrix<double, 12, 12> M = Eigen::Matrix<double, 12, 12>::Zero();
        double L = length;
        double m = rho * area * L;

        // Lump half translational mass to each node (u_x, u_y, u_z)
        for (int d = 0; d < 3; ++d) {
            M(d, d) = m * 0.5;
            M(d + 6, d + 6) = m * 0.5;
        }

        // Approximate rotational inertia (very approximate, for dynamics only)
        // Use a heuristic Irot = m * L^2 / 12 (like rod's inertia) split half to each node
        double Irot = m * L * L / 12.0;
        for (int d = 3; d < 6; ++d) {
            M(d, d) = Irot * 0.5;
            M(d + 6, d + 6) = Irot * 0.5;
        }

        return M;
    }


//class CantileverBeam3D {
//private:
//    double length;           // Total beam length (m)
//    double E;                // Young's modulus (Pa)
//    double nu;               // Poisson's ratio
//    double rho;              // Density (kg/m^3)
//    double area;             // Cross-sectional area (m^2)
//    double Iyy, Izz, J;      // section properties
//    int numElements;         // Number of finite elements
//    double baseOutd, IzzBase;
//    Eigen::Vector3d endPointLoad; // Point load at the free end (N), 3D vector
//
//    std::vector<BeamElement3D> elements;
//    Eigen::MatrixXd globalStiffnessMatrix;
//    Eigen::MatrixXd globalMassMatrix;
//    Eigen::MatrixXd globalDampingMatrix;
//    Eigen::VectorXd forceVectorMag;
//    Eigen::VectorXd forceVector;
//    Eigen::VectorXd u; // displacement
//    Eigen::VectorXd v; // velocity
//
//
//public:
    CantileverBeam3D::CantileverBeam3D(double _length, double _E, double _nu, double _rho, double _area,
        int _numElements, const Eigen::Vector3d& _endPointLoad,
        double _outDia, double _inDia, double _taper)
        : length(_length), E(_E), nu(_nu), rho(_rho), area(_area),
        numElements(_numElements), endPointLoad(_endPointLoad) {

        double elementLength = length / numElements;
        baseOutd = _outDia;

        // Build elements with tapered properties along the span
        for (int i = 0; i < numElements; ++i) {
            double taper_factor = 1.0 - _taper * double(i) / double(numElements);
            double outD = _outDia * taper_factor;
            double inD = _inDia * taper_factor;
            double A = calculateHollowTubeArea(outD, inD);
            double I = calculateHollowTubeMomentOfInertia(outD, inD); // Iyy == Izz for circular
            if (i == 0) {
                IzzBase = I;
            }
            double Jp = 2.0 * I; // polar moment for circular tube Jp = Iyy + Izz = 2*I
            double G = E / (2.0 * (1.0 + nu));
            elements.emplace_back(elementLength, E, G, I, I, Jp, rho, A);
        }


        // Global DOFs: 6 per node
        int totalDOFs = 6 * (numElements + 1);
        globalStiffnessMatrix = Eigen::MatrixXd::Zero(totalDOFs, totalDOFs);
        globalMassMatrix = Eigen::MatrixXd::Zero(totalDOFs, totalDOFs);
        forceVectorMag = Eigen::VectorXd::Zero(totalDOFs);
        forceVector = forceVectorMag;

        assembleGlobalMatrices();
        applyBoundaryConditions();
        applyLoads();
    }

    void CantileverBeam3D::assembleGlobalMatrices() {
        for (int e = 0; e < numElements; ++e) {
            Eigen::Matrix<double, 12, 12> localK = elements[e].getStiffnessMatrix();
            Eigen::Matrix<double, 12, 12> localM = elements[e].getMassMatrix();

            int startDOF = 6 * e;
            for (int r = 0; r < 12; ++r) {
                for (int c = 0; c < 12; ++c) {
                    int gr = startDOF + r;
                    int gc = startDOF + c;
                    globalStiffnessMatrix(gr, gc) += localK(r, c);
                    globalMassMatrix(gr, gc) += localM(r, c);
                }
            }
        }
    }

    void CantileverBeam3D::applyBoundaryConditions() {
        // Fix all 6 DOFs at the first node (cantilever root) using heavy penalty
        double penalty = 1e15;
        for (int i = 0; i < 6; ++i) {
            globalStiffnessMatrix(i, i) = penalty;
            // zero out coupling to be safe
            for (int j = 0; j < (int)globalStiffnessMatrix.rows(); ++j) {
                if (j != i) {
                    globalStiffnessMatrix(i, j) = 0.0;
                    globalStiffnessMatrix(j, i) = 0.0;
                }
            }
        }
    }

    void CantileverBeam3D::applyLoads() {
        // Apply the end point load at the last node translations (u_x,u_y,u_z)
        int lastNode = numElements;
        int base = 6 * lastNode;
        forceVector(base + 0) = endPointLoad(0);
        forceVector(base + 1) = endPointLoad(1);
        forceVector(base + 2) = endPointLoad(2);
        forceVectorMag = forceVector;
    }




    // Sleep for a given number of milliseconds (calls std::this_thread::sleep_for)
    static inline void sleep_ms(std::uint64_t ms) {
        std::this_thread::sleep_for(std::chrono::milliseconds(ms));
    }

    // Monotonic time in milliseconds (use for elapsed-time measurements)
    static inline std::uint64_t now_millis_steady() noexcept {
        using namespace std::chrono;
        return static_cast<std::uint64_t>(
            duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count()
            );
    }

    void sync(double dt)
    {
        static uint64_t last_run = 0;
        const uint64_t interval_ms = (uint64_t)(dt * 1000.0);
        uint64_t now = now_millis_steady();
        if (now - last_run < interval_ms) {
            uint64_t naptime_ms = interval_ms - (now - last_run);
            sleep_ms(naptime_ms);
        }
    }

    //Just draw the beam
    void CantileverBeam3D::draw(openGLframe& graphics)
    {
        // visualize current state: reconstruct full displacement vector
        int totalDOFs = 6 * (numElements + 1);
        int activeDOFs = totalDOFs - 6; // excluding fixed first node
        Eigen::VectorXd fullDisp = Eigen::VectorXd::Zero(totalDOFs);
        fullDisp.tail(activeDOFs) = u;

        std::vector<LineSegment> lineSegments;
        float sc = 1.8f / float(numElements + 1);
        for (int n = 0; n < numElements; ++n) {
            float x1 = sc * n - 0.9f + static_cast<float>(fullDisp(6 * n));
            float x2 = sc * (n + 1) - 0.9f + static_cast<float>(fullDisp(6 * (n + 1)));
            //float x1 = 1.0f * static_cast<float>(fullDisp(6 * n ));
            //float x2 = 1.0f * static_cast<float>(fullDisp(6 * (n + 1)));
            float y1 = 1.0f * static_cast<float>(fullDisp(6 * n + 1));
            float y2 = 1.0f * static_cast<float>(fullDisp(6 * (n + 1) + 1));
            double drot_y = fullDisp(6 * (n + 1) + 4) - fullDisp(6 * n + 4);
            double drot_z = fullDisp(6 * (n + 1) + 5) - fullDisp(6 * n + 5);
            double bend = static_cast<float>(std::sqrt(drot_y * drot_y + drot_z * drot_z) * 1000.0);
            double stress = abs(static_cast<float>(fullDisp(6 * n)) - static_cast<float>(fullDisp(6 * (n + 1)))) * 500.0;
            RGBi color = valueToHeatmapColor(bend + stress);
            lineSegments.emplace_back(x1, y1, x2, y2, color.r / 255.0f, color.g / 255.0f, color.b / 255.0f);
        }
        //graphics.CLS(RGBi(0, 0, 0));
        graphics.drawLines(lineSegments);
        //graphics.swap();
    }

    void CantileverBeam3D::showOnScreen(openGLframe& graphics, double dt)
    {
        if (dt > 0) sync(dt);
        graphics.CLS(RGBi(0, 0, 0));
        draw(graphics);
        graphics.swap();
    }

    //void CantileverBeam3D::visualize(openGLframe& graphics, Eigen::VectorXd& u, int _numElements, double dt)
    //{
    //    if (dt > 0) sync(dt);
    //    // visualize current state: reconstruct full displacement vector
    //    int totalDOFs = 6 * (_numElements + 1);
    //    int activeDOFs = totalDOFs - 6; // excluding fixed first node
    //    Eigen::VectorXd fullDisp = Eigen::VectorXd::Zero(totalDOFs);
    //    fullDisp.tail(activeDOFs) = u;

    //    std::vector<LineSegment> lineSegments;
    //    float sc = 1.8f / float(numElements + 1);
    //    for (int n = 0; n < numElements; ++n) {
    //        float x1 = sc * n - 0.9f + static_cast<float>(fullDisp(6 * n));
    //        float x2 = sc * (n + 1) - 0.9f + static_cast<float>(fullDisp(6 * (n + 1)));
    //        //float x1 = 1.0f * static_cast<float>(fullDisp(6 * n ));
    //        //float x2 = 1.0f * static_cast<float>(fullDisp(6 * (n + 1)));
    //        float y1 = 1.0f * static_cast<float>(fullDisp(6 * n + 1));
    //        float y2 = 1.0f * static_cast<float>(fullDisp(6 * (n + 1) + 1));
    //        double drot_y = fullDisp(6 * (n + 1) + 4) - fullDisp(6 * n + 4);
    //        double drot_z = fullDisp(6 * (n + 1) + 5) - fullDisp(6 * n + 5);
    //        double bend = static_cast<float>(std::sqrt(drot_y * drot_y + drot_z * drot_z) * 1000.0);
    //        double stress = abs(static_cast<float>(fullDisp(6 * n)) - static_cast<float>(fullDisp(6 * (n + 1)))) * 500.0;
    //        RGBi color = valueToHeatmapColor(bend + stress);
    //        lineSegments.emplace_back(x1, y1, x2, y2, color.r / 255.0f, color.g / 255.0f, color.b / 255.0f);
    //    }
    //    graphics.CLS(RGBi(0, 0, 0));
    //    graphics.drawLines(lineSegments);
    //    graphics.swap();
    //}


    //Option 4: Measure Curvature Directly
    //    Instead of displacement, measure the curvature or strain at the base :
    // M = EI⋅κ = EI⋅d2vdx2
    // M = EI \cdot \kappa = EI \cdot \frac{ d ^ 2v }{dx ^ 2}
    //    Then :

    // Fy = dMdx
    // F_y = \frac{ dM }{dx}

    void CantileverBeam3D::solveStaticDisplacement() {
        int totalDOFs = 6 * (numElements + 1);
        int activeDOFs = totalDOFs - 6; // excluding fixed first node

        Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Eigen::VectorXd Factive = forceVector.tail(activeDOFs);

        //Eigen::VectorXd displacements = Kactive.colPivHouseholderQr().solve(Factive);
        u = Kactive.colPivHouseholderQr().solve(Factive);

        Eigen::VectorXd fullDisp = Eigen::VectorXd::Zero(totalDOFs);
        fullDisp.tail(activeDOFs) = u;

        std::cout << "Static Displacement Analysis Results (3D):" << std::endl;
        std::cout << "------------------------------------------" << std::endl;

        std::cout << "Forces transferred at base:" << std::endl;
        int mp = 1;//measurement point
        double unitLength = mp * length / numElements;
        double Fx = (fullDisp(6 * mp * (1)) / unitLength) * E * area;
        std::cout << "Fx:" << Fx << std::endl;
        double uy1 = fullDisp(6 * mp + 1);
        double th_y1 = fullDisp(6 * mp + 4);
        double th_z1 = fullDisp(6 * mp + 5);
        double uy2 = fullDisp(6 * mp * 2 + 1);
        double uy3 = fullDisp(6 * mp * 3 + 1);
        double ddy1 = uy1 / unitLength;
        double ddy2 = (uy2 - uy1) / unitLength;
        double ddy3 = (uy3 - uy2) / unitLength;
        double ddysqr1 = (ddy2 - ddy1) / unitLength;
        double ddysqr2 = (ddy3 - ddy2) / unitLength;
        double mt = E * IzzBase * ddysqr1;
        double Fy2 = E * IzzBase * (ddysqr2 - ddysqr1) / unitLength;
        double Ftest = E * IzzBase * (uy3 - 3 * uy2 + 3 * uy1) / pow(unitLength, 3);
        //double Fy = uy1 * 12 * E * IzzBase / pow(length, 3) + 6*mp * E * IzzBase * th_z1 / pow(length, 2);
        double Fy = 6 * E * IzzBase / pow(unitLength, 3) * (uy1 - unitLength * th_z1 / 2);
        double m = E * IzzBase * (uy2 - 2 * uy1) / pow(unitLength, 2);
        std::cout << "Fy:" << Fy << "(Izz=" << IzzBase << ")" << std::endl;
        std::cout << "Fy2:" << Fy2 << std::endl;
        double Fz = fullDisp(6 * mp + 2) * 12 * E * IzzBase / pow(unitLength, 3) + 6 * E * IzzBase * fullDisp(6 * mp + 4) / pow(unitLength, 2);
        std::cout << "Fz:" << Fz << std::endl;
        double Mz = 2 * E * IzzBase / pow(unitLength, 2) * (2 * th_z1 - 3 * uy1 / unitLength);
        std::cout << "Mt:" << mt << std::endl;
        double Mz2 = E * IzzBase * th_z1 / unitLength; // curvature
        double My = E * IzzBase * th_y1 / unitLength; //curvature
        std::cout << "(method 2) My:" << My << std::endl;
        std::cout << "(method 2) Mz:" << Mz2 << std::endl;
        std::cout << std::endl;

        for (int n = 0; n <= numElements; ++n) {
            double pos = double(n) * (length / numElements);
            std::cout << "Node " << n << " (x = " << pos << " m):" << std::endl;
            std::cout << "  ux: " << fullDisp(6 * n + 0) << " m" << std::endl;
            if (n < numElements)
                std::cout << "  delta ux: " << fullDisp(6 * (n + 1) + 0) - fullDisp(6 * n + 0) << " m" << std::endl;

            std::cout << "  uy: " << fullDisp(6 * n + 1) << " m" << std::endl;
            std::cout << "  uz: " << fullDisp(6 * n + 2) << " m" << std::endl;
            std::cout << "  rot_x: " << fullDisp(6 * n + 3) << " rad" << std::endl;
            std::cout << "  rot_y: " << fullDisp(6 * n + 4) << " rad" << std::endl;
            std::cout << "  rot_z: " << fullDisp(6 * n + 5) << " rad" << std::endl;
            std::cout << std::endl;
        }


    }

    void CantileverBeam3D::visualize(openGLframe& graphics)
    {
        // Simple visualization: plot uy vs x
        showOnScreen(graphics, -1);
    }

    void CantileverBeam3D::solveFrequencyAnalysis(int numModes) {
        int totalDOFs = 6 * (numElements + 1);
        int activeDOFs = totalDOFs - 6;
        Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Eigen::MatrixXd Mactive = globalMassMatrix.bottomRightCorner(activeDOFs, activeDOFs);

        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(Kactive, Mactive);

        std::cout << "Natural Frequencies (3D):" << std::endl;
        for (int i = 0; i < std::min(numModes, activeDOFs); ++i) {
            double omega2 = solver.eigenvalues()(i);
            if (omega2 <= 0) continue;
            double freq = std::sqrt(omega2) / (2.0 * M_PI);
            std::cout << "Mode " << (i + 1) << ": " << freq << " Hz" << std::endl;
        }
    }


    int totalDOFs, activeDOFs;
    Eigen::MatrixXd Kactive, Mactive, Factive, Cactive;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    Eigen::VectorXd acc;
    Eigen::LLT<Eigen::MatrixXd> lltOfLHS;
    double timeStep;


    void CantileverBeam3D::setupTimeDomainSimulation(double _timeStep, double dampingRatio)
    {
        totalDOFs = 6 * (numElements + 1);
        activeDOFs = totalDOFs - 6;
        Kactive = globalStiffnessMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Mactive = globalMassMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Factive = forceVector.tail(activeDOFs);
        timeStep = _timeStep;

        // Rayleigh damping using two lowest modes (if available)
        //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        solver.compute(Kactive, Mactive);
        double omega1 = 0.0, omega2 = 0.0;
        if (solver.eigenvalues().size() >= 2) {
            omega1 = std::sqrt(std::max(0.0, solver.eigenvalues()(0)));
            omega2 = std::sqrt(std::max(0.0, solver.eigenvalues()(1)));
        }
        else if (solver.eigenvalues().size() == 1) {
            omega1 = std::sqrt(std::max(0.0, solver.eigenvalues()(0)));
            omega2 = 2.0 * omega1;
        }
        else {
            omega1 = 1.0; omega2 = 2.0;
        }
        double alpha = dampingRatio * 2.0 * omega1 * omega2 / (omega1 + omega2);
        double beta = dampingRatio * 2.0 / (omega1 + omega2);
        Cactive = alpha * Mactive + beta * Kactive;

        // Newmark parameters
        const double gamma = 0.5;
        const double beta_nb = 0.25;
        u = Eigen::VectorXd::Zero(activeDOFs); // displacement
        v = Eigen::VectorXd::Zero(activeDOFs); // velocity
        acc = Mactive.colPivHouseholderQr().solve(Factive - Kactive * u - Cactive * v);

        Eigen::MatrixXd LHS = Mactive + gamma * timeStep * Cactive + beta_nb * timeStep * timeStep * Kactive;
        //Eigen::LLT<Eigen::MatrixXd> lltOfLHS(LHS);
        lltOfLHS = Eigen::LLT<Eigen::MatrixXd>(LHS);
    }

    void CantileverBeam3D::stepForward(double timeStep)
    {
        static int step = 0;
        double time = step++ * timeStep;
        const double gamma = 0.5;
        const double beta_nb = 0.25;

        // --- Update time-varying forces here ---
        double amp = 1.0;
        if (time < 10)
            amp = sin(10 * time);
        else if (time < 20)
            amp = 1.0;
        else
            amp = 0;

        int lastNode = numElements;
        int base = 6 * lastNode;
        forceVector(base + 0) = amp * endPointLoad(0);
        forceVector(base + 1) = amp * endPointLoad(1);
        forceVector(base + 2) = amp * endPointLoad(2);

        Factive = forceVector.tail(activeDOFs);

        // visualize current state: reconstruct full displacement vector
        //visualize(graphics, u, numElements, timeStep);

        // Newmark predictor
        Eigen::VectorXd uPred = u + timeStep * v + timeStep * timeStep * (0.5 - beta_nb) * acc;
        Eigen::VectorXd vPred = v + timeStep * (1.0 - gamma) * acc;

        Eigen::VectorXd RHS = Factive - Kactive * uPred - Cactive * vPred;

        Eigen::VectorXd deltaA = lltOfLHS.solve(RHS);

        acc = deltaA;
        v = vPred + gamma * timeStep * acc;
        u = uPred + beta_nb * timeStep * timeStep * acc;
    }

    void CantileverBeam3D::simulateTimeDomain(openGLframe& graphics, double duration, double _timeStep, double _dampingRatio)
    {
        setupTimeDomainSimulation(_timeStep, _dampingRatio);
        int numSteps = static_cast<int>(duration / _timeStep);
        for (int step = 0; step <= numSteps; ++step) {
            stepForward(_timeStep);
            showOnScreen(graphics, _timeStep);
        }
    }


    // Time integration (Newmark) simplified for 3D DOFs; uses lumped/approx mass so it's stable-ish.
    void CantileverBeam3D::simulateTimeDomain2(openGLframe& graphics, double duration, double timeStep, double dampingRatio) {
        int totalDOFs = 6 * (numElements + 1);
        int activeDOFs = totalDOFs - 6;
        Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Eigen::MatrixXd Mactive = globalMassMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Eigen::VectorXd Factive = forceVector.tail(activeDOFs);

        // Rayleigh damping using two lowest modes (if available)
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        solver.compute(Kactive, Mactive);
        double omega1 = 0.0, omega2 = 0.0;
        if (solver.eigenvalues().size() >= 2) {
            omega1 = std::sqrt(std::max(0.0, solver.eigenvalues()(0)));
            omega2 = std::sqrt(std::max(0.0, solver.eigenvalues()(1)));
        }
        else if (solver.eigenvalues().size() == 1) {
            omega1 = std::sqrt(std::max(0.0, solver.eigenvalues()(0)));
            omega2 = 2.0 * omega1;
        }
        else {
            omega1 = 1.0; omega2 = 2.0;
        }
        double alpha = dampingRatio * 2.0 * omega1 * omega2 / (omega1 + omega2);
        double beta = dampingRatio * 2.0 / (omega1 + omega2);
        Eigen::MatrixXd Cactive = alpha * Mactive + beta * Kactive;

        // Newmark parameters
        double gamma = 0.5;
        double beta_nb = 0.25;
        Eigen::VectorXd u = Eigen::VectorXd::Zero(activeDOFs); // displacement
        Eigen::VectorXd v = Eigen::VectorXd::Zero(activeDOFs); // velocity
        Eigen::VectorXd acc = Mactive.colPivHouseholderQr().solve(Factive - Kactive * u - Cactive * v);

        Eigen::MatrixXd LHS = Mactive + gamma * timeStep * Cactive + beta_nb * timeStep * timeStep * Kactive;
        Eigen::LLT<Eigen::MatrixXd> lltOfLHS(LHS);

        int numSteps = static_cast<int>(duration / timeStep);

        for (int step = 0; step <= numSteps; ++step) {
            double time = step * timeStep;

            // --- Update time-varying forces here ---
            double amp = 1.0;
            if (time < 10)
                amp = sin(10 * time);
            else if (time < 20)
                amp = 1.0;
            else
                amp = 0;

            int lastNode = numElements;
            int base = 6 * lastNode;
            forceVector(base + 0) = amp * endPointLoad(0);
            forceVector(base + 1) = amp * endPointLoad(1);
            forceVector(base + 2) = amp * endPointLoad(2);

            Factive = forceVector.tail(activeDOFs);

            // visualize current state: reconstruct full displacement vector
            showOnScreen(graphics, timeStep);

            // Newmark predictor
            Eigen::VectorXd uPred = u + timeStep * v + timeStep * timeStep * (0.5 - beta_nb) * acc;
            Eigen::VectorXd vPred = v + timeStep * (1.0 - gamma) * acc;

            Eigen::VectorXd RHS = Factive - Kactive * uPred - Cactive * vPred;

            Eigen::VectorXd deltaA = lltOfLHS.solve(RHS);

            acc = deltaA;
            v = vPred + gamma * timeStep * acc;
            u = uPred + beta_nb * timeStep * timeStep * acc;
            //u[0] = sin(time);
            //std::cout << u[1] << std::endl;
        }
    }





CantileverBeam3D make_beam()
{

    // Hollow aluminum tube parameters
    double length = 45.0 * 12.0 / 39.37;          // Length (ft -> m)
    double outerDiameter = 8.0 / 39.37;           // 8" -> m
    double innerDiameter = outerDiameter - 2.0 * thickness;
    double E = 69e9;                              // Young's modulus for aluminum (Pa)
    double nu = 0.33;                             // Poisson's ratio
    double rho = 2700.0;                          // Density (kg/m^3)
    double taper = .5;// 75. / 100.;

    double area;
    double I;
    double J; // polar moment for circular tube

    area = calculateHollowTubeArea(outerDiameter, innerDiameter);
    I = calculateHollowTubeMomentOfInertia(outerDiameter, innerDiameter);
    J = 2.0 * I; // polar moment for circular tube

    int numElements = 20;
    // Point load at free end in global (Fx, Fy, Fz). Original used vertical N; here we apply in Y direction
    Eigen::Vector3d pointLoad(0, 100.0, 0.0); // convert lbf to N and apply in Y

    std::cout << "3D Finite Element Cantilever (Beam) - Hollow Aluminum Tube" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << "Beam length: " << length << " m" << std::endl;
    std::cout << "Base Outer diameter: " << outerDiameter << " m" << std::endl;
    std::cout << "Base Inner diameter: " << innerDiameter << " m" << std::endl;
    std::cout << "End Outer diameter: " << outerDiameter * (1 - taper) << " m" << std::endl;
    std::cout << "End Inner diameter: " << innerDiameter * (1 - taper) << " m" << std::endl;
    std::cout << "Base Cross-sectional area: " << area << " m^2" << std::endl;
    std::cout << "Moment of inertia (Iyy=Izz): " << I << " m^4" << std::endl;
    std::cout << "Polar moment approx J: " << J << " m^4" << std::endl;
    std::cout << "Young's modulus: " << E << " Pa" << std::endl;
    std::cout << "Poisson ratio: " << nu << std::endl;
    std::cout << "Density: " << rho << " kg/m^3" << std::endl;
    std::cout << "Point load at free end (N): (" << pointLoad.transpose() << ")" << std::endl;
    std::cout << "Number of elements: " << numElements << std::endl << std::endl;

    CantileverBeam3D beam(length, E, nu, rho, area, numElements, pointLoad, outerDiameter, innerDiameter, taper);
    return beam;
}

#ifndef BEAM2_SUBMODULE
int main()
{
    openGLframe graphics;
    graphics.setupGL();

#else
int beam2_init(openGLframe &graphics)
#endif
{
    //// Hollow aluminum tube parameters
    //double length = 45.0 * 12.0 / 39.37;          // Length (ft -> m)
    //double outerDiameter = 8.0 / 39.37;           // 8" -> m
    //double innerDiameter = outerDiameter - 2.0 * thickness;
    //double E = 69e9;                              // Young's modulus for aluminum (Pa)
    //double nu = 0.33;                             // Poisson's ratio
    //double rho = 2700.0;                          // Density (kg/m^3)
    //double taper = .5;// 75. / 100.;

    //double area;
    //double I;
    //double J; // polar moment for circular tube

    //area = calculateHollowTubeArea(outerDiameter, innerDiameter);
    //I = calculateHollowTubeMomentOfInertia(outerDiameter, innerDiameter);
    //J = 2.0 * I; // polar moment for circular tube

    //int numElements = 20;
    //// Point load at free end in global (Fx, Fy, Fz). Original used vertical N; here we apply in Y direction
    //Eigen::Vector3d pointLoad(  0, 100.0, 0.0); // convert lbf to N and apply in Y

    //std::cout << "3D Finite Element Cantilever (Beam) - Hollow Aluminum Tube" << std::endl;
    //std::cout << "==========================================================" << std::endl;
    //std::cout << "Beam length: " << length << " m" << std::endl;
    //std::cout << "Base Outer diameter: " << outerDiameter << " m" << std::endl;
    //std::cout << "Base Inner diameter: " << innerDiameter << " m" << std::endl;
    //std::cout << "End Outer diameter: " << outerDiameter*(1-taper) << " m" << std::endl;
    //std::cout << "End Inner diameter: " << innerDiameter*(1-taper) << " m" << std::endl;
    //std::cout << "Base Cross-sectional area: " << area << " m^2" << std::endl;
    //std::cout << "Moment of inertia (Iyy=Izz): " << I << " m^4" << std::endl;
    //std::cout << "Polar moment approx J: " << J << " m^4" << std::endl;
    //std::cout << "Young's modulus: " << E << " Pa" << std::endl;
    //std::cout << "Poisson ratio: " << nu << std::endl;
    //std::cout << "Density: " << rho << " kg/m^3" << std::endl;
    //std::cout << "Point load at free end (N): (" << pointLoad.transpose() << ")" << std::endl;
    //std::cout << "Number of elements: " << numElements << std::endl << std::endl;

    //CantileverBeam3D beam(length, E, nu, rho, area, numElements, pointLoad, outerDiameter, innerDiameter, taper);

    CantileverBeam3D beam = make_beam();

    // Static analysis
    //beam.solveStaticDisplacement();
    //beam.showOnScreen(graphics);

    std::cout << std::endl;

    // Modal analysis (first 6 modes)
    //beam.solveFrequencyAnalysis(6);

    std::cout << std::endl;

    // Time domain simulation
    beam.simulateTimeDomain(graphics, 60.0, 1/30.0);

    graphics.waitForCompletion();
    graphics.closeGL();

    return 0;
}