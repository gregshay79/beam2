#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <cstdint>
#include <thread>
#include <chrono>
#include "beam2.h"

#ifndef BEAM2_SUBMODULE
//#include "../opengl3.h"
#include "../heatmap.h"
#else
//#include "../opengl3.h"
#include "../heatmap.h"
#endif

#define M_PI 3.141592653589793238462643383


double thickness = 0.525 / 39.37;              // 1/4" -> m
//extern double mast_height;


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
        K(0, 0) = EA_L; 
        K(0, 6) = -EA_L;
        K(6, 0) = -EA_L; 
        K(6, 6) = EA_L;

        // Torsion
        double GJ_L = G * J / L;
        K(3, 3) = GJ_L; 
        K(3, 9) = -GJ_L;
        K(9, 3) = -GJ_L; 
        K(9, 9) = GJ_L;

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

    CantileverBeam3D::CantileverBeam3D(double _length, double _E, double _nu, double _rho, double _area,
        int _numElements, const Eigen::Vector3d& _endPointLoad,
        double _outDia, double _inDia, double _taper)
        : length(_length), E(_E), nu(_nu), rho(_rho), area(_area),
        numElements(_numElements), endPointLoad(_endPointLoad) {

        // Global DOFs: 6 per node
        totalDOFs = 6 * (numElements + 1);
        activeDOFs = totalDOFs;// -6;
        elementLength = length / numElements;
        baseOutd = _outDia;

        x = Eigen::VectorXd::Zero(activeDOFs); // displacement
        v = Eigen::VectorXd::Zero(activeDOFs); // velocity
        acc = Eigen::VectorXd::Zero(activeDOFs); // velocity

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
            Eigen::Vector3d undeformed_pos(i * elementLength, 0.0, 0.0);
            x.segment<3>(i * 6) = undeformed_pos;

        }

        globalStiffnessMatrix = Eigen::MatrixXd::Zero(totalDOFs, totalDOFs);
        globalMassMatrix = Eigen::MatrixXd::Zero(totalDOFs, totalDOFs);

        assembleGlobalMatrices();
        
        K_coupling = globalStiffnessMatrix.block(totalDOFs - activeDOFs, 0, activeDOFs, 6);
        M_coupling = globalMassMatrix.block(totalDOFs-activeDOFs, 0, activeDOFs, 6);
        //K_base_to_active_original = globalStiffnessMatrix.block(0, 6, 6, activeDOFs);
        //M_base_to_active_original = globalMassMatrix.block(0, 6, 6, activeDOFs);
        //K_base_base_original = globalStiffnessMatrix.block(0, 0, 6, 6);
        //M_base_base_original = globalMassMatrix.block(0, 0, 6, 6);
        
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

    Eigen::VectorXd CantileverBeam3D::applyEndpointLoad(Eigen::Vector3d endPointLoad) {
        // Apply the end point load at the last node translations (u_x,u_y,u_z)
        int lastNode = numElements;
        int base = 6 * lastNode;
        Eigen::VectorXd _forceVector = Eigen::VectorXd::Zero(totalDOFs);
        _forceVector(base + 0) = endPointLoad(0);
        _forceVector(base + 1) = endPointLoad(1);
        _forceVector(base + 2) = endPointLoad(2);
        return _forceVector;
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

    //wait until at least 'dt' time since the last sync call.
    //If time is later than last+dt, go ahead immediately.
    void sync(double dt)
    {
        static uint64_t last_run = 0;
        const uint64_t interval_ms = (uint64_t)(dt * 1000.0);
        uint64_t now = now_millis_steady();
        if (now - last_run < interval_ms) {
            uint64_t naptime_ms = interval_ms - (now - last_run);
            sleep_ms(naptime_ms);
        }
        last_run = now;
    }

    //void CantileverBeam3D::setBaseOffset(Eigen::Vector3d _offset)
    //{
    //    u(3+0) =  _offset(0);
    //    u(3+1) =  _offset(1);
    //    u(3+2) =  _offset(2);
    //    u(6 + 1) = elementLength * sin(_offset(2));
    //}
    // 
    //Just draw the beam
    void CantileverBeam3D::draw(openGLframe& graphics)
    {
        // visualize current state: reconstruct full displacement vector
        Eigen::VectorXd fullDisp = Eigen::VectorXd::Zero(totalDOFs);
        fullDisp.tail(activeDOFs) = x; // u is of length ActiveDOF

        std::vector<LineSegment> lineSegments;
        //float sc = (float)mast_height / float(numElements + 1);
        float sc = (float)length / float(numElements + 1);
        for (int n = 0; n < numElements; ++n) {
            //float y1 = (float)yb + (sc * n + static_cast<float>(fullDisp(6 * n)));
            //float y2 = (float)yb + (sc * (n + 1) + static_cast<float>(fullDisp(6 * (n + 1))));
            float y1 = (static_cast<float>(fullDisp(6 * n)));
            float y2 = (static_cast<float>(fullDisp(6 * (n + 1))));

            float x1 = (1.0f * static_cast<float>(fullDisp(6 * n + 1))) ;
            float x2 = (1.0f * static_cast<float>(fullDisp(6 * (n + 1) + 1)));

            x1 = (float)(x1 * gxscale + goffset);
            x2 = (float)(x2 * gxscale + goffset);
            y1 = (float)(y1 * gyscale + goffset);
            y2 = (float)(y2 * gyscale + goffset);

            double drot_y = fullDisp(6 * (n + 1) + 4) - fullDisp(6 * n + 4);
            double drot_z = fullDisp(6 * (n + 1) + 5) - fullDisp(6 * n + 5);
            double bend = static_cast<float>(std::sqrt(drot_y * drot_y + drot_z * drot_z) * 1000.0);
            double stress = abs(static_cast<float>(fullDisp(6 * n)) - static_cast<float>(fullDisp(6 * (n + 1)))) * 500.0;
            RGBi color = valueToHeatmapColor((bend /* + stress*/)/3.0 );
            lineSegments.emplace_back(x1, y1, x2, y2, color.r / 255.0f, color.g / 255.0f, color.b / 255.0f);
        }
        //graphics.CLS(RGBi(0, 0, 0));
        graphics.drawLines(lineSegments);
        //graphics.swap();
    }

    //show on screen.  If dt is given, synchronize screen refresh with dt.
    void CantileverBeam3D::showOnScreen(openGLframe& graphics, double dt)
    {
        if (dt > 0) sync(dt);
        graphics.CLS(RGBi{ 0, 0, 0 });
        draw(graphics);
        graphics.swap();
    }

    //Option 4: Measure Curvature Directly
    //    Instead of displacement, measure the curvature or strain at the base :
    // M = EI⋅κ = EI⋅d2vdx2
    // M = EI \cdot \kappa = EI \cdot \frac{ d ^ 2v }{dx ^ 2}
    //    Then :

    // Fy = dMdx
    // F_y = \frac{ dM }{dx}

    void CantileverBeam3D::solveStaticDisplacement(Eigen::VectorXd& forceVector) {
        Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Eigen::VectorXd Factive = forceVector.tail(activeDOFs);

        //Eigen::VectorXd displacements = Kactive.colPivHouseholderQr().solve(Factive);
        Eigen::VectorXd u = Eigen::VectorXd::Zero(activeDOFs);
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
        double mt = E * IzzBase * ddysqr1; //this is most accurate measure of momentum at base
        double Fy2 = E * IzzBase * (ddysqr2 - ddysqr1) / unitLength; //This is most accurate measure of force
        //double Ftest = E * IzzBase * (uy3 - 3 * uy2 + 3 * uy1) / pow(unitLength, 3);
        //double Fy = uy1 * 12 * E * IzzBase / pow(length, 3) + 6*mp * E * IzzBase * th_z1 / pow(length, 2);
        //double Fy = 6 * E * IzzBase / pow(unitLength, 3) * (uy1 - unitLength * th_z1 / 2);
        //double m = E * IzzBase * (uy2 - 2 * uy1) / pow(unitLength, 2);
        //std::cout << "Fy:" << Fy << "(Izz=" << IzzBase << ")" << std::endl;
        std::cout << "Fy2:" << Fy2 << std::endl;
        double Fz = fullDisp(6 * mp + 2) * 12 * E * IzzBase / pow(unitLength, 3) + 6 * E * IzzBase * fullDisp(6 * mp + 4) / pow(unitLength, 2);
        std::cout << "Fz:" << Fz << std::endl;
        double Mz = 2 * E * IzzBase / pow(unitLength, 2) * (2 * th_z1 - 3 * uy1 / unitLength);
        std::cout << "Mt:" << mt << std::endl;
        double Mz2 = E * IzzBase * th_z1 / unitLength; // curvature method.  Close, but somewhat low.
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
        showOnScreen(graphics); //unsynchronized
    }

    void CantileverBeam3D::solveFrequencyAnalysis(int numModes) {
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

    void CantileverBeam3D::setupTimeDomainSimulation(double _timeStep, double dampingRatio)
    {
        // Extract active submatrices
        Kactive = globalStiffnessMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Mactive = globalMassMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        
        // Use original coupling matrices (saved before boundary conditions modified the matrix)
        //K_coupling = K_coupling_original;
        //M_coupling = M_coupling_original;
        //K_base_to_active = K_base_to_active_original;
        //M_base_to_active = M_base_to_active_original;
        //K_base_base = K_base_base_original;
        //M_base_base = M_base_base_original;
        
        // Initialize base state to zero
        //u_base = Eigen::VectorXd::Zero(6);
        //v_base = Eigen::VectorXd::Zero(6);
        //acc_base = Eigen::VectorXd::Zero(6);
        
        Eigen::VectorXd forceVector = Eigen::VectorXd::Zero(activeDOFs);
        Factive = forceVector.tail(activeDOFs);
        timeStep = _timeStep;

        // Rayleigh damping using two lowest modes (if available)
        //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;
        //solver.compute(Kactive, Mactive);
        double omega1 = 0.0, omega2 = 0.0;
        //if (solver.eigenvalues().size() >= 2) {
        //    omega1 = std::sqrt(std::max(0.0, solver.eigenvalues()(0)));
        //    omega2 = std::sqrt(std::max(0.0, solver.eigenvalues()(1)));
        //}
        //else if (solver.eigenvalues().size() == 1) {
        //    omega1 = std::sqrt(std::max(0.0, solver.eigenvalues()(0)));
        //    omega2 = 2.0 * omega1;
        //}
        //else {
            omega1 = 1.0; omega2 = 2.0;
        //}
        double alpha = dampingRatio * 2.0 * omega1 * omega2 / (omega1 + omega2);
        double beta = dampingRatio * 2.0 / (omega1 + omega2);
        Cactive = alpha * Mactive + beta * Kactive;
        
        // Damping coupling matrices
        C_coupling = alpha * M_coupling + beta * K_coupling;
        //C_base_to_active = alpha * M_base_to_active + beta * K_base_to_active;

        // Newmark parameters
        const double gamma = 0.5;
        const double beta_nb = 0.25;
        //Note: positions x already computed in constructor
        //u = Eigen::VectorXd::Zero(activeDOFs); // displacement
        v = Eigen::VectorXd::Zero(activeDOFs); // velocity
        //acc = Mactive.colPivHouseholderQr().solve(Factive - Kactive * u - Cactive * v);
        acc = Eigen::VectorXd::Zero(activeDOFs);
        
        Eigen::MatrixXd LHS = Mactive + gamma * timeStep * Cactive + beta_nb * timeStep * timeStep * Kactive;
        lltOfLHS = Eigen::LLT<Eigen::MatrixXd>(LHS);
    }

/* test modifying forces 
    // --- Update time-varying forces here ---
    static int step = 0;
    double time = step++ * timeStep;
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

*/
// Helper: Convert small-angle rotations (rx, ry, rz) to rotation matrix
// For large angles, use proper Euler angle conversion
Eigen::Matrix3d rotationMatrixFromAngles(const Eigen::Vector3d& angles) {
    double rx = angles(0);
    double ry = angles(1);
    double rz = angles(2);
    
    // ZYX Euler angles (commonly used in beam theory)
    Eigen::Matrix3d Rx, Ry, Rz;
    Rx << 1, 0, 0,
          0, cos(rx), -sin(rx),
          0, sin(rx), cos(rx);
    
    Ry << cos(ry), 0, sin(ry),
          0, 1, 0,
          -sin(ry), 0, cos(ry);
    
    Rz << cos(rz), -sin(rz), 0,
          sin(rz), cos(rz), 0,
          0, 0, 1;
    
    return Rz * Ry * Rx;  // Order matters for large angles
}

// Extract angles from rotation matrix (inverse operation)
Eigen::Vector3d anglesFromRotationMatrix(const Eigen::Matrix3d& R) {
    double ry = asin(-R(2,0));
    double rx = atan2(R(2,1), R(2,2));
    double rz = atan2(R(1,0), R(0,0));
    return Eigen::Vector3d(rx, ry, rz);
}

void CantileverBeam3D::stepForward(double timeStep, Eigen::VectorXd& forceVector)
{
    const double gamma = 0.5;
    const double beta_nb = 0.25;

    // Local storage for transformed variables: Transform all nodes to body-fixed frame
    Eigen::VectorXd u_body = Eigen::VectorXd::Zero(activeDOFs);
    Eigen::VectorXd v_body = v;
    Eigen::VectorXd acc_body = acc;
    Eigen::VectorXd forceVector_body = forceVector;  // ADDED: force vector in body frame

    // Convert to body-fixed (relative) coordinates
    Eigen::Vector3d base_pos = x.head(3);
    Eigen::Vector3d base_rot = x.segment<3>(3);
    // Get base velocities and accelerations (from first node)
    Eigen::Vector3d base_vel_linear = v.segment<3>(0);
    Eigen::Vector3d base_vel_angular = v.segment<3>(3);
    Eigen::Vector3d base_acc_linear = acc.segment<3>(0);
    Eigen::Vector3d base_acc_angular = acc.segment<3>(3);

    Eigen::Matrix3d R_base = rotationMatrixFromAngles(base_rot);
    Eigen::Matrix3d R_base_inv = R_base.transpose();  // Inverse rotation
    
    // Loop through all nodes and transform to body-fixed frame
    for (int i = 0; i < numElements + 1; ++i) {
        int idx = 6 * i;
        
        // Transform positions to body frame, then subtract undeformed position
        Eigen::Vector3d pos_global = x.segment<3>(idx);
        Eigen::Vector3d pos_relative = pos_global - base_pos;
        Eigen::Vector3d pos_in_body_frame = R_base_inv * pos_relative;
        
        // Undeformed position in body frame (beam runs along x-axis)
        Eigen::Vector3d undeformed_pos(i * elementLength, 0.0, 0.0);
        u_body.segment<3>(idx) = pos_in_body_frame - undeformed_pos;
        
        // Transform rotations (small angle approximation: simple subtraction)
        Eigen::Vector3d rot_global = x.segment<3>(idx + 3);
        u_body.segment<3>(idx + 3) = rot_global - base_rot;
        
        // Transform velocities - make relative to base THEN rotate
        Eigen::Vector3d vel_relative = v.segment<3>(idx) - base_vel_linear;
        v_body.segment<3>(idx) = R_base_inv * vel_relative;
        
        Eigen::Vector3d angvel_relative = v.segment<3>(idx + 3) - base_vel_angular;
        v_body.segment<3>(idx + 3) = R_base_inv * angvel_relative;
        
        // Transform accelerations - make relative to base THEN rotate
        Eigen::Vector3d acc_relative = acc.segment<3>(idx) - base_acc_linear;
        acc_body.segment<3>(idx) = R_base_inv * acc_relative;
        
        Eigen::Vector3d angacc_relative = acc.segment<3>(idx + 3) - base_acc_angular;
        acc_body.segment<3>(idx + 3) = R_base_inv * angacc_relative;
        
        // Transform forces and moments to body frame
        Eigen::Vector3d force_global = forceVector.segment<3>(idx);
        forceVector_body.segment<3>(idx) = R_base_inv * force_global;
        
        Eigen::Vector3d moment_global = forceVector.segment<3>(idx + 3);
        forceVector_body.segment<3>(idx + 3) = R_base_inv * moment_global;
    }
    
    // Set base node to origin in body-fixed frame
    if ((u_body.head(6).norm() > .01) || (v_body.head(6).norm() > .01) || (acc_body.head(6).norm() > 0.01)) {
        printf("problem!\n");
    }
    //u_body.head(6).setZero(); //should already be 0
    //v_body.head(6).setZero(); // should already be 0
    //acc_body.head(6).setZero();//should already be 0

    // Use transformed forces
    Factive = forceVector_body.tail(activeDOFs);

    // Newmark integration in body-fixed frame
    Eigen::VectorXd uPred = u_body + timeStep * v_body + timeStep * timeStep * (0.5 - beta_nb) * acc_body;
    Eigen::VectorXd vPred = v_body + timeStep * (1.0 - gamma) * acc_body;

    Eigen::VectorXd RHS = Factive - Kactive * uPred - Cactive * vPred;
    Eigen::VectorXd deltaA = lltOfLHS.solve(RHS);

    acc_body = deltaA;
    v_body = vPred + gamma * timeStep * acc_body;
    u_body = uPred + beta_nb * timeStep * timeStep * acc_body;

    // Transform back to global coordinates
    for (int i = 0; i < numElements + 1; ++i) {
        int idx = 6 * i;
        
        // Transform positions back: add undeformed position, rotate, then add base position
        Eigen::Vector3d deformation_body = u_body.segment<3>(idx);
        Eigen::Vector3d undeformed_pos(i * elementLength, 0.0, 0.0);
        Eigen::Vector3d pos_in_body_frame = deformation_body + undeformed_pos;
        x.segment<3>(idx) = R_base * pos_in_body_frame + base_pos;
        
        // Transform rotations back (small angle approximation: simple addition)
        Eigen::Vector3d rot_body = u_body.segment<3>(idx + 3);
        x.segment<3>(idx + 3) = rot_body + base_rot;
        
        // Transform velocities back - rotate THEN add base velocity
        Eigen::Vector3d vel_body = v_body.segment<3>(idx);
        v.segment<3>(idx) = R_base * vel_body + base_vel_linear;
        
        Eigen::Vector3d angvel_body = v_body.segment<3>(idx + 3);
        v.segment<3>(idx + 3) = R_base * angvel_body + base_vel_angular;
        
        // Transform accelerations back - rotate THEN add base acceleration
        Eigen::Vector3d acc_body_linear = acc_body.segment<3>(idx);
        acc.segment<3>(idx) = R_base * acc_body_linear + base_acc_linear;
        
        Eigen::Vector3d acc_body_angular = acc_body.segment<3>(idx + 3);
        acc.segment<3>(idx + 3) = R_base * acc_body_angular + base_acc_angular;
    }

    
}

    // void CantileverBeam3D::stepForward(double timeStep, Eigen::VectorXd& forceVector)
    // {
    //     const double gamma = 0.5;
    //     const double beta_nb = 0.25;

    //     //convert to relative coordinates
    //     base_root = u.head(6);
    //     for (int i = 0; i < totalDOFs; i += 6) {
    //         u.block(i,0,6, 1) -= base_root;
    //     }

    //     Factive = forceVector.tail(activeDOFs);

    //     // Newmark predictor
    //     Eigen::VectorXd uPred = u + timeStep * v + timeStep * timeStep * (0.5 - beta_nb) * acc;
    //     Eigen::VectorXd vPred = v + timeStep * (1.0 - gamma) * acc;

    //     Eigen::VectorXd RHS = Factive - Kactive * uPred - Cactive * vPred;

    //     Eigen::VectorXd deltaA = lltOfLHS.solve(RHS);

    //     acc = deltaA;
    //     v = vPred + gamma * timeStep * acc;
    //     u = uPred + beta_nb * timeStep * timeStep * acc;

    //     //convert back to global coordinates
    //     for (int i = 0; i < totalDOFs; i += 6) {
    //         u.block(i,0,6,1) += base_root;
    //     }
    // }


    //void CantileverBeam3D::stepForwardWithBaseMotion(double timeStep, 
    //                                                  const Eigen::VectorXd& externalForces)
    //{
    //    const double gamma = 0.5;
    //    const double beta_nb = 0.25;
    //    
    //    // Compute forces on active DOFs due to base motion
    //    Eigen::VectorXd F_base_motion = -K_coupling * u_base 
    //                                   - C_coupling * v_base 
    //                                   - M_coupling * acc_base;

    //    //double fx = F_base_motion(0);
    //    //double fy = F_base_motion(1);
    //    //double fz = F_base_motion(2);
    //    
    //    // Extract active DOF forces from external forces (externalForces is totalDOFs size)
    //    // Base DOFs (0-5) are ignored since base motion is prescribed separately
    //    Eigen::VectorXd F_external_active = externalForces.tail(activeDOFs);
    //    
    //    
    //    // Combine external forces with base motion forces
    //    Eigen::VectorXd F_total = F_external_active + F_base_motion ;;
    //    
    //    Factive = F_total;
    //    
    //    // Newmark integration (same as before)
    //    Eigen::VectorXd uPred = u + timeStep * v + timeStep * timeStep * (0.5 - beta_nb) * acc;
    //    Eigen::VectorXd vPred = v + timeStep * (1.0 - gamma) * acc;
    //    
    //    Eigen::VectorXd RHS = Factive - Kactive * uPred - Cactive * vPred;
    //    Eigen::VectorXd deltaA = lltOfLHS.solve(RHS);
    //    
    //    acc = deltaA;
    //    v = vPred + gamma * timeStep * acc;
    //    u = uPred + beta_nb * timeStep * timeStep * acc;
    //}

    //Eigen::VectorXd CantileverBeam3D::getBaseReactionForces() const
    //{
    //    // Reaction forces = forces needed to maintain base motion
    //    // F_reaction = K(base, base) * u_base + K(base, active) * u_active
    //    //            + C(base, base) * v_base + C(base, active) * v_active  
    //    //            + M(base, base) * acc_base + M(base, active) * acc_active
    //    
    //    Eigen::VectorXd F_reaction = Eigen::VectorXd::Zero(6);
    //    
    //    // Forces from base motion itself
    //    //F_reaction += K_base_base * u_base + M_base_base * acc_base;
    //    
    //    // Forces from coupling to active DOFs
    //    F_reaction += K_base_to_active * u;
    //    F_reaction += M_base_to_active * acc;
    //    
    //    // Add damping coupling
    //    F_reaction += C_base_to_active * v;
    //    
    //    return F_reaction;
    //}

    void CantileverBeam3D::simulateTimeDomain(openGLframe& graphics, double duration, double _timeStep, double _dampingRatio)
    {
        setupTimeDomainSimulation(_timeStep, _dampingRatio);
        Eigen::VectorXd forceVector = applyEndpointLoad(endPointLoad);
        int numSteps = static_cast<int>(duration / _timeStep);
        for (int step = 0; step <= numSteps; ++step) {
            double time = step * timeStep;
            double amp = 1.0;
            if (time < 20)
                amp = sin(10 * time);
            else if (time < 40)
                amp = 1.0;
            else
                amp = 0;

            amp *= 10;
            int lastNode = numElements;
            int base = 6 * lastNode;
            forceVector(base + 0) = amp * endPointLoad(0);
            forceVector(base + 1) = amp * endPointLoad(1);
            forceVector(base + 2) = amp * endPointLoad(2);

            //add a weak spring to hold the base at the origin
            Eigen::Vector3d bs = x.head(3);
            Eigen::Vector3d force = -1*bs.normalized();

            double dist = bs.norm();
            force *= 1000.0 * dist; //weak spring force

            // add a damping force proportional to the square of the velocity
            Eigen::Vector3d vec = v.head(3);
            double dampingForceMag = -10.0 * vec.squaredNorm();
            force += dampingForceMag * vec.normalized();

            forceVector.head(3) = force;

            //Add a moment holding the base to angle 0
            Eigen::Vector3d angle_vec = x.segment<3>(3);
            double angle_mag = angle_vec.norm();
            angle_vec.normalize();
            Eigen::Vector3d righting_moment = -10000.0 * angle_mag * angle_vec;
            forceVector.segment<3>(3) = righting_moment;


            stepForward(_timeStep, forceVector);
            showOnScreen(graphics, _timeStep);
        }
    }


    // Time integration (Newmark) simplified for 3D DOFs; uses lumped/approx mass so it's stable-ish.
    void CantileverBeam3D::simulateTimeDomain2(openGLframe& graphics, double duration, double timeStep, double dampingRatio) {
        Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Eigen::MatrixXd Mactive = globalMassMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Eigen::VectorXd forceVector = applyEndpointLoad(endPointLoad);
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
            //modify forcevector
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
    gxscale *= 0.8;
    gyscale *= 0.8;
    goffset = 0.0;
    graphics.setupGL();

#else
int beam2_init(openGLframe &graphics)
{
#endif

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
    Eigen::Vector3d endPointLoad(  0, 40000.0, 0.0); // convert lbf to N and apply in Y

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
    Eigen::VectorXd forceVector = beam.applyEndpointLoad(endPointLoad);

    //beam.solveStaticDisplacement(forceVector);
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