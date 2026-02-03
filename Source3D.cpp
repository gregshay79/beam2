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

    void CantileverBeam3D::AttachBase(double mass, Eigen::Vector3d rotational_inertia) {

        // Add large mass to base node (node 0)
        for (int i = 0; i < 3; ++i) {  // translational DOFs (x, y, z)
            globalMassMatrix(i, i) += mass;
        }
        
        // Add large rotational inertia to base node
        for (int i = 3; i < 6; ++i) {  // rotational DOFs (rx, ry, rz)
            globalMassMatrix(i, i) += rotational_inertia(i-3);
        }

        //Recompute coupling matrices
        //K_coupling = globalStiffnessMatrix.block(totalDOFs - activeDOFs, 0, activeDOFs, 6);
        //M_coupling = globalMassMatrix.block(totalDOFs-activeDOFs, 0, activeDOFs, 6);
    }

    CantileverBeam3D::CantileverBeam3D(double _length, double _E, double _nu, double _rho, double _area,
        int _numElements, double _outDia, double _inDia, double _taper)
        : length(_length), E(_E), nu(_nu), rho(_rho), area(_area),
        numElements(_numElements) {

        // Global DOFs: 6 per node
        totalDOFs = 6 * (numElements + 1);
        activeDOFs = totalDOFs;// -6;
        elementLength = length / numElements;
        baseOutd = _outDia;

        x = Eigen::VectorXd::Zero(activeDOFs); // position
        ref_pos = Eigen::VectorXd::Zero(activeDOFs); // reference position
        v = Eigen::VectorXd::Zero(activeDOFs); // velocity
        acc = Eigen::VectorXd::Zero(activeDOFs); // velocity
        origin_displacement = Eigen::Vector3d::Zero();
        //origin_orientation = Eigen::Vector3d::Zero();

        // Build elements with tapered properties along the span
        int i;
        for (i = 0; i < numElements; ++i) {
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
            ref_pos.segment<3>(i * 6) = undeformed_pos;
        }
        Eigen::Vector3d undeformed_pos(i * elementLength, 0.0, 0.0); //compute last node
        ref_pos.segment<3>(i * 6) = undeformed_pos;

        globalStiffnessMatrix = Eigen::MatrixXd::Zero(totalDOFs, totalDOFs);
        globalMassMatrix = Eigen::MatrixXd::Zero(totalDOFs, totalDOFs);

        assembleGlobalMatrices();
        
        //K_coupling = globalStiffnessMatrix.block(totalDOFs - activeDOFs, 0, activeDOFs, 6);
        //M_coupling = globalMassMatrix.block(totalDOFs-activeDOFs, 0, activeDOFs, 6);

        double angle = 0;// M_PI / 16;
        for (i = 0; i < numElements + 1; ++i) {
            x.segment<3>(i * 6) = Eigen::Vector3d(sin(angle) * i * elementLength, cos(angle) * i * elementLength, 0);
            x.segment<3>(i * 6 + 3) = Eigen::Vector3d(0, 0, 0);
        }
        
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

    // Relocate the position vector x to be at the new position.
    void CantileverBeam3D::jam_position(Eigen::Vector3d new_position)
    {
        Eigen::Vector3d base = x.segment<3>(0);
        for (int n = 0; n < numElements+1; ++n) {
            x.segment<3>(6*n) -= base;
            x.segment<3>(6*n) += new_position;
            //x(6 * n + 4) = M_PI / 32;
        }
    }

    //Just draw the beam
    void CantileverBeam3D::draw(openGLframe& graphics)
    {
        // visualize current state: reconstruct full displacement vector
        Eigen::VectorXd fullDisp = Eigen::VectorXd::Zero(totalDOFs);
        fullDisp.tail(activeDOFs) = x; // u is of length ActiveDOF

        std::vector<LineSegment> lineSegments;

        for (int n = 0; n < numElements; ++n) {
            float x1 = (static_cast<float>(fullDisp(6 * n)));
            float x2 = (static_cast<float>(fullDisp(6 * (n + 1))));

            float y1 = (static_cast<float>(fullDisp(6 * n + 1))) ;
            float y2 = (static_cast<float>(fullDisp(6 * (n + 1) + 1)));

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
        graphics.drawLines(lineSegments);
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
        // base0_global_equilibrium = { 0,0,0 };
        // base1_global_equilibrium = { 0.0, elementLength, 0.0 };
        Eigen::Vector3d rot_inertia = { 8000,18000,18000 };//in beam relative space
        AttachBase(7500, rot_inertia);

        // Extract active submatrices AFTER adding springs
        Kactive = globalStiffnessMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Mactive = globalMassMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        
        Eigen::VectorXd forceVector = Eigen::VectorXd::Zero(activeDOFs);
        Factive = forceVector.tail(activeDOFs);
        timeStep = _timeStep;
        origin_displacement = { 0,0,0 };
        //origin_orientation = { 0,0,0 };

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
        //C_coupling = alpha * M_coupling + beta * K_coupling;

        // Newmark parameters
        gamma = 0.5;
        beta_nb = 0.25;
        //Note: positions x already computed in constructor
        v = Eigen::VectorXd::Zero(activeDOFs); // velocity
        //acc = Mactive.colPivHouseholderQr().solve(Factive - Kactive * u - Cactive * v);
        acc = Eigen::VectorXd::Zero(activeDOFs);
        
        Eigen::MatrixXd LHS = Mactive + gamma * timeStep * Cactive + beta_nb * timeStep * timeStep * Kactive;
        lltOfLHS = Eigen::LLT<Eigen::MatrixXd>(LHS);
    }

// Compute rotation matrix that rotates x-axis to align with given direction
Eigen::Matrix3d rotationFromXAxisTo(const Eigen::Vector3d& direction) {
    Eigen::Vector3d dir = direction.normalized();
    Eigen::Vector3d x_axis(1, 0, 0);
    
    // If direction is already along x-axis, return identity
    if ((dir - x_axis).norm() < 1e-6) {
        return Eigen::Matrix3d::Identity();
    }
    
    // If direction is opposite x-axis, rotate 180° around z
    if ((dir + x_axis).norm() < 1e-6) {
        Eigen::Matrix3d R;
        R << -1, 0, 0,
              0, -1, 0,
              0, 0, 1;
        return R;
    }
    
    // General case: use Rodrigues' rotation formula
    Eigen::Vector3d v = x_axis.cross(dir);  // rotation axis
    double s = v.norm();                     // sin(angle)
    double c = x_axis.dot(dir);              // cos(angle)
    
    Eigen::Matrix3d vx;  // skew-symmetric matrix of v
    vx << 0, -v(2), v(1),
          v(2), 0, -v(0),
          -v(1), v(0), 0;
    
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + vx + vx * vx * ((1 - c) / (s * s));
    
    return R;
}

// Compute rotation matrix that rotates x-axis to align with given direction
Eigen::Matrix3d rotationFromYAxisTo(const Eigen::Vector3d& direction) {
    Eigen::Vector3d dir = direction.normalized();
    Eigen::Vector3d y_axis(0, 1, 0);

    // If direction is already along x-axis, return identity
    if ((dir - y_axis).norm() < 1e-6) {
        return Eigen::Matrix3d::Identity();
    }

    // If direction is opposite x-axis, rotate 180° around z
    if ((dir + y_axis).norm() < 1e-6) {
        Eigen::Matrix3d R;
        R << -1, 0, 0,
            0, -1, 0,
            0, 0, 1;
        return R;
    }

    // General case: use Rodrigues' rotation formula
    Eigen::Vector3d v = y_axis.cross(dir);  // rotation axis
    double s = v.norm();                     // sin(angle)
    double c = y_axis.dot(dir);              // cos(angle)

    Eigen::Matrix3d vx;  // skew-symmetric matrix of v
    vx << 0, -v(2), v(1),
        v(2), 0, -v(0),
        -v(1), v(0), 0;

    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + vx + vx * vx * ((1 - c) / (s * s));

    return R;
}

// Extract Euler angles from rotation matrix (ZYX convention: yaw-pitch-roll)
// Avoids gimbal lock by detecting singularities
Eigen::Vector3d eulerAnglesZYX(const Eigen::Matrix3d& R)
{
    Eigen::Vector3d angles;
    
    // Check for gimbal lock (pitch near ±90°)
    double sin_pitch = -R(2, 0);
    
    const double GIMBAL_LOCK_THRESHOLD = 0.99999;
    
    if (std::abs(sin_pitch) > GIMBAL_LOCK_THRESHOLD)
    {
        // Gimbal lock case: pitch is near ±90°
        // In this case, yaw and roll are not uniquely defined
        // We set yaw to 0 and compute roll as the combination
        angles(0) = 0.0;  // yaw (rotation about Z)
        angles(1) = std::asin(std::clamp(sin_pitch, -1.0, 1.0));  // pitch (rotation about Y)
        
        if (sin_pitch > 0)  // pitch = +90°
        {
            angles(2) = std::atan2(-R(0, 1), R(1, 1));  // roll (rotation about X)
        }
        else  // pitch = -90°
        {
            angles(2) = std::atan2(R(0, 1), R(1, 1));
        }
    }
    else
    {
        // Normal case: no gimbal lock
        angles(0) = std::atan2(R(1, 0), R(0, 0));  // yaw (rotation about Z)
        angles(1) = std::asin(std::clamp(sin_pitch, -1.0, 1.0));  // pitch (rotation about Y)
        angles(2) = std::atan2(R(2, 1), R(2, 2));  // roll (rotation about X)
    }
    
    return angles;  // Returns [yaw, pitch, roll] in radians
}

Eigen::VectorXd CantileverBeam3D::getOrientation()
{

    return x.segment<6>(0);
    Eigen::VectorXd result = Eigen::VectorXd::Zero(6);
    //// Determine beam's current orientation from the first element position and direction.
    //result.segment<3>(0) = x.segment<3>(0);

    ////Wrong: compute the angles relative to the axis, using the first two nodes of the beam.
    ////Eigen::Vector3d node0_pos = x.segment<3>(0);
    ////Eigen::Vector3d node1_pos = x.segment<3>(6);
    ////Eigen::Vector3d beam_direction = (node1_pos - node0_pos).normalized();
    ////Eigen::Matrix3d R_beam = rotationFromYAxisTo(beam_direction);
    ////Eigen::Vector3d angles = eulerAnglesZYX(R_beam);
    ////result.segment<3>(3) = angles;

    //result.segment<3>(3) = x.segment<3>(3);

    //return result;
}

Eigen::VectorXd CantileverBeam3D::getBaseVelocity()
{
    return v.head(6);
}

// Convert XYZ Euler angles to rotation matrix
// Rotation order: first X, then Y, then Z
Eigen::Matrix3d rotationMatrixFromAnglesXYZ(const Eigen::Vector3d& angles) {
    double angle_x = angles(0);  // rotation about X-axis (roll)
    double angle_y = angles(1);  // rotation about Y-axis (pitch)
    double angle_z = angles(2);  // rotation about Z-axis (yaw)
    
    // Individual rotation matrices
    Eigen::Matrix3d Rx, Ry, Rz;
    
    // Rotation about X-axis
    Rx << 1, 0,            0,
          0, cos(angle_x), -sin(angle_x),
          0, sin(angle_x),  cos(angle_x);
    
    // Rotation about Y-axis
    Ry << cos(angle_y),  0, sin(angle_y),
          0,             1, 0,
          -sin(angle_y), 0, cos(angle_y);
    
    // Rotation about Z-axis
    Rz << cos(angle_z), -sin(angle_z), 0,
          sin(angle_z),  cos(angle_z), 0,
          0,             0,            1;
    
    // XYZ order: apply X first, then Y, then Z
    // R = Rz * Ry * Rx
    return Rz * Ry * Rx;
}

void CantileverBeam3D::stepForward(double timeStep, Eigen::VectorXd& forceVector)
{
    // Determine beam's current orientation from first element direction
    // Eigen::Vector3d node0_pos = x.segment<3>(0);
    // Eigen::Vector3d node1_pos = x.segment<3>(6);
    // Eigen::Vector3d beam_direction = (node1_pos - node0_pos).normalized();

    //Eigen::Vector3d base_pos = x.segment<3>(0);
    Eigen::Vector3d base_pos = x.segment<3>(0) - origin_displacement; // subtract the displacement of the base. WOW! this was correct!

    //No... Find the rotation from the x frame of reference to the beam frame of reference.
    //This method below loses the y axis yaw
    // Eigen::Matrix3d R_beam = rotationFromXAxisTo(beam_direction);
    // Eigen::Matrix3d R_inv = R_beam.transpose();

    // Compute rotation matrix from the angles in the first element,
    // to transform from the world frame to the beam frame (rotation 0,0,0 is the beam frame).
    Eigen::Vector3d base_angles = x.segment<3>(3);
 
    Eigen::Matrix3d R_beam = rotationMatrixFromAnglesXYZ(base_angles);

    // Create π/2 rotation about Z-axis
    Eigen::Matrix3d R_z_90;
    R_z_90 << 0, -1, 0,
            1,  0, 0,
            0,  0, 1;

    Eigen::Matrix3d R_z_minus90;
    R_z_minus90 <<  0,  1, 0,
                    -1,  0, 0,
                    0,  0, 1;
    // Apply transpose of beam rotation first, then Z rotation in beam's local frame
    R_beam = R_z_minus90 * R_beam.transpose();
    Eigen::Matrix3d R_inv = R_beam;
        
    // Transform positions, velocities, and forces to beam-aligned frame
    Eigen::VectorXd x_rot(activeDOFs);
    Eigen::VectorXd v_rot = v;
    Eigen::VectorXd acc_rot = acc;
    Eigen::VectorXd f_rot = forceVector;
    
    for (int i = 0; i < numElements + 1; ++i) {
        int idx = 6 * i;
        // Rotate positions to beam frame
        x_rot.segment<3>(idx) = R_inv * (x.segment<3>(idx)-base_pos);
        
        // Transform rotations to beam-aligned frame (treat angles as vectors for small angles)
        x_rot.segment<3>(idx + 3) = R_inv * x.segment<3>(idx + 3);
        
        // Rotate velocities and accelerations
        v_rot.segment<3>(idx) = R_inv * v.segment<3>(idx);
        acc_rot.segment<3>(idx) = R_inv * acc.segment<3>(idx);
        
        // Rotate forces and moments
        f_rot.segment<3>(idx) = R_inv * forceVector.segment<3>(idx);
        f_rot.segment<3>(idx + 3) = R_inv * forceVector.segment<3>(idx + 3);
    }
    
    // NOW compute displacement in beam-aligned frame
    Eigen::VectorXd u_rot = x_rot - ref_pos;
    
    Factive = f_rot.tail(activeDOFs);
        
    // Newmark integration in beam-aligned frame (where stiffness matrix is correct)
    Eigen::VectorXd uPred = u_rot + timeStep * v_rot + timeStep * timeStep * (0.5 - beta_nb) * acc_rot;
    Eigen::VectorXd vPred = v_rot + timeStep * (1.0 - gamma) * acc_rot;
    
    Eigen::VectorXd RHS = Factive - Kactive * uPred - Cactive * vPred;
    Eigen::VectorXd deltaA = lltOfLHS.solve(RHS);
    
    acc_rot = deltaA;
    v_rot = vPred + gamma * timeStep * acc_rot;
    u_rot = uPred + beta_nb * timeStep * timeStep * acc_rot;
    
    // Compute positions in beam frame from displacements
    Eigen::VectorXd x_rot_new = ref_pos + u_rot;


    // update the base angles after the new force deformation.
    // The new force deformation angles are small, so can directly add.
    base_angles += u_rot.segment<3>(3);

    // Rotate back to world space.
    // First rotate by -90 degrees about Z-axis, then rotate by the updated base angles.
    Eigen::Matrix3d R_beam_updated = rotationMatrixFromAnglesXYZ(base_angles);
    R_beam_updated = R_beam_updated * R_z_90;

    // Find the beam origin displacement in global space
    origin_displacement = R_beam_updated * u_rot.segment<3>(0);
    
    // Transform back to global frame
    for (int i = 0; i < numElements + 1; ++i) {
        int idx = 6 * i;
        // Rotate positions back
        x.segment<3>(idx) = R_beam_updated * (x_rot_new.segment<3>(idx)) + base_pos;
        
        // Transform rotations back to global frame
        x.segment<3>(idx + 3) = R_beam_updated * x_rot_new.segment<3>(idx + 3);
        
        // Rotate velocities and accelerations back
        v.segment<3>(idx) = R_beam_updated * v_rot.segment<3>(idx);
        acc.segment<3>(idx) = R_beam_updated * acc_rot.segment<3>(idx);
    }

    // static int count = 0;
    // count++;
    // if (count >= 30) {
    //     count = 0;
    //     // Compute spring reaction forces
    //     //Eigen::Vector3d base_equilibrium(0.0, 0.0, 0.0);
    //     //Eigen::Vector3d node1_equilibrium(0.0, elementLength, 0.0);

    //     Eigen::Vector3d force_on_holding_spring = k_holding * (x.head(3)-base0_global_equilibrium);
    //     Eigen::Vector3d force_on_righting_spring = k_righting * (x.segment<3>(6) - base1_global_equilibrium);

    //     // Print or store these forces
    //     std::cout << std::endl << "Force on base spring: ["
    //             << force_on_holding_spring(0) << ", "
    //             << force_on_holding_spring(1) << ", "
    //             << force_on_holding_spring(2) << "]" << std::endl;
                
    //     std::cout << "Force on node1 spring: [" 
    //             << force_on_righting_spring(0) << ", "
    //             << force_on_righting_spring(1) << ", "
    //             << force_on_righting_spring(2) << "]" << std::endl;    
    // }
}
 
    void CantileverBeam3D::simulateTimeDomain(openGLframe& graphics, double duration, double _timeStep, double _dampingRatio)
    {
        setupTimeDomainSimulation(_timeStep, _dampingRatio);
        Eigen::Vector3d pointLoad = { 400,0,0 }; //100 lbs force
        Eigen::VectorXd forceVector = Eigen::VectorXd::Zero(activeDOFs);
        int numSteps = static_cast<int>(duration / _timeStep);
        for (int step = 0; step <= numSteps; ++step) {
            double time = step * timeStep;
            double amp = 1.0;
            forceVector = Eigen::VectorXd::Zero(activeDOFs);
            if (time < 10) {
                amp = 5.0;
                if (time < 5)
                    amp = sin(2 * time);
            }
            else {
                amp = 0;
            }

            
            //amp = 0;

            if (time > 20) amp = 0;
            
            int pushedNode = numElements;
            int tip = 6 * pushedNode;
            forceVector.segment<3>(tip) += amp * pointLoad;


            // Springs are now implicit (in stiffness matrix)
            // Only add: (1) equilibrium position forces, (2) velocity damping

            // Explicit velocity damping (can't be implicit)
            // Base velocity damping
            //double base_damping = -300.0;
            //forceVector.segment<3>(0) += base_damping * v.head(3);
            
            // Swing damping at midpoint
            //int midNode = numElements / 2;
            //double swing_damping = -300.0;
            //forceVector.segment<3>(6 * midNode) += swing_damping * v.segment<3>(6 * midNode);

            stepForward(_timeStep, forceVector);
            showOnScreen(graphics, _timeStep);
        }
    }


    // Time integration (Newmark) simplified for 3D DOFs; uses lumped/approx mass so it's stable-ish.
    void CantileverBeam3D::simulateTimeDomain2(openGLframe& graphics, double duration, double timeStep, double dampingRatio) {
        Eigen::MatrixXd Kactive = globalStiffnessMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Eigen::MatrixXd Mactive = globalMassMatrix.bottomRightCorner(activeDOFs, activeDOFs);
        Eigen::VectorXd forceVector = Eigen::VectorXd::Zero(activeDOFs);
        Eigen::VectorXd Factive = forceVector.tail(activeDOFs);
        Eigen::Vector3d pointLoad = { 10000,0,0 };

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
            forceVector.segment<3>(base) = amp * pointLoad;

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
    //Eigen::Vector3d pointLoad(0, 100.0, 0.0); // convert lbf to N and apply in Y

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
    //std::cout << "Point load at free end (N): (" << pointLoad.transpose() << ")" << std::endl;
    std::cout << "Number of elements: " << numElements << std::endl << std::endl;

    CantileverBeam3D beam(length, E, nu, rho, area, numElements, outerDiameter, innerDiameter, taper);
    return beam;
}

#ifndef BEAM2_SUBMODULE
int main()
{
    openGLframe graphics;
    gxscale *= 0.8;
    gyscale *= 0.8;
    goffset = -0.25;
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
    //Eigen::Vector3d endPointLoad(0,0.0, 0.0); // convert lbf to N and apply in Y

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
    //Eigen::VectorXd forceVector = beam.applyEndpointLoad(endPointLoad);
    Eigen::VectorXd forceVector = Eigen::VectorXd::Zero(beam.DOF());
    Eigen::Vector3d pointLoad = { 10000,0,0 };
    forceVector.segment<3>(beam.DOF() - 6) = pointLoad;
    //beam.solveStaticDisplacement(forceVector);
    //beam.showOnScreen(graphics);

    std::cout << std::endl;

    // Modal analysis (first 6 modes)
    //beam.solveFrequencyAnalysis(6);

    std::cout << std::endl;

    // Time domain simulation
    beam.simulateTimeDomain(graphics, 300.0, 1/30.0);

    graphics.waitForCompletion();
    graphics.closeGL();

    return 0;
}