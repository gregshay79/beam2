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
        Eigen::Vector3d undeformed_pos(i * elementLength, 0.0, 0.0); //insert last node
        ref_pos.segment<3>(i * 6) = undeformed_pos;

        //Now set initial position of nodes of beam
        double angle = 5.0*M_PI/180.0;
        for (i = 0; i < numElements + 1; ++i) {
//            x.segment<3>(i * 6) = Eigen::Vector3d(0.0, i * elementLength, 0);
//            x.segment<3>(i * 6 + 3) = Eigen::Vector3d(0.0, 0.0, M_PI/2);
            x.segment<3>(i * 6) = Eigen::Vector3d(sin(angle)*i * elementLength,cos(angle)*i * elementLength, 0);
            x.segment<3>(i * 6 + 3) = Eigen::Vector3d(0.0, 0.0, M_PI/4);
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

    //Eigen::VectorXd CantileverBeam3D::applyEndpointLoad(Eigen::Vector3d endPointLoad) {
    //    // Apply the end point load at the last node translations (u_x,u_y,u_z)
    //    int lastNode = numElements;
    //    int base = 6 * lastNode;
    //    Eigen::VectorXd _forceVector = Eigen::VectorXd::Zero(totalDOFs);
    //    _forceVector(base + 0) = endPointLoad(0);
    //    _forceVector(base + 1) = endPointLoad(1);
    //    _forceVector(base + 2) = endPointLoad(2);
    //    return _forceVector;
    //}




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


    //Just draw the beam
    void CantileverBeam3D::draw(openGLframe& graphics)
    {
        // visualize current state: reconstruct full displacement vector
        Eigen::VectorXd fullDisp = Eigen::VectorXd::Zero(totalDOFs);
        fullDisp.tail(activeDOFs) = x; // u is of length ActiveDOF

        std::vector<LineSegment> lineSegments;
        //float sc = (float)mast_height / float(numElements + 1);
        //float sc = (float)length / float(numElements + 1);
        for (int n = 0; n < numElements; ++n) {
            //float y1 = (float)yb + (sc * n + static_cast<float>(fullDisp(6 * n)));
            //float y2 = (float)yb + (sc * (n + 1) + static_cast<float>(fullDisp(6 * (n + 1))));
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
        gamma = 0.5;
        beta_nb = 0.25;
        //Note: positions x already computed in constructor
        //u = Eigen::VectorXd::Zero(activeDOFs); // displacement
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

void CantileverBeam3D::stepForward(double timeStep, Eigen::VectorXd& forceVector)
{
    // Determine beam's current orientation from first element direction
    Eigen::Vector3d node0_pos = x.segment<3>(0);
    Eigen::Vector3d node1_pos = x.segment<3>(6);
    Eigen::Vector3d beam_direction = (node1_pos - node0_pos).normalized();
    
    // Compute rotation from reference (x-axis) to current beam direction
    Eigen::Matrix3d R_beam = rotationFromXAxisTo(beam_direction);
    Eigen::Matrix3d R_inv = R_beam.transpose();
    
    // Transform positions, velocities, and forces to beam-aligned frame
    Eigen::VectorXd x_rot(activeDOFs);
    Eigen::VectorXd v_rot = v;
    Eigen::VectorXd acc_rot = acc;
    Eigen::VectorXd f_rot = forceVector;
    
    for (int i = 0; i < numElements + 1; ++i) {
        int idx = 6 * i;
        // Rotate positions to beam frame
        x_rot.segment<3>(idx) = R_inv * x.segment<3>(idx);
        
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
    
    // Transform back to global frame
    for (int i = 0; i < numElements + 1; ++i) {
        int idx = 6 * i;
        // Rotate positions back
        x.segment<3>(idx) = R_beam * x_rot_new.segment<3>(idx);
        
        // Transform rotations back to global frame
        x.segment<3>(idx + 3) = R_beam * x_rot_new.segment<3>(idx + 3);
        
        // Rotate velocities and accelerations back
        v.segment<3>(idx) = R_beam * v_rot.segment<3>(idx);
        acc.segment<3>(idx) = R_beam * acc_rot.segment<3>(idx);
    }
}

    // void CantileverBeam3D::stepForward(double timeStep, Eigen::VectorXd& forceVector)
    // {


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
        Eigen::Vector3d pointLoad = { 100,0,0 };
        Eigen::VectorXd forceVector = Eigen::VectorXd::Zero(activeDOFs);
        int numSteps = static_cast<int>(duration / _timeStep);
        for (int step = 0; step <= numSteps; ++step) {
            double time = step * timeStep;
            double amp = 1.0;
            if (time < 40) {
                amp = 5.0;
                if (time < 20)
                    amp = sin(2 * time);
            }
            else {
                amp = 0;
            }

            amp *= 20;
            int pushedNode = numElements / 2;
            int tip = 6 * pushedNode;
            forceVector.segment<3>(tip) = amp * pointLoad;
            //forceVector(tip + 0) = amp * pointLoad(0);
            //forceVector(tip + 1) = amp * pointLoad(1);
            //forceVector(tip + 2) = amp * pointLoad(2);

            //add a weak spring to hold the base at the origin
            Eigen::Vector3d bs = x.head(3);
            Eigen::Vector3d force = -1 * bs.normalized();

            double dist = bs.norm();
            force *= 50000.0 * dist; //base spring force

            // add a damping force proportional to the velocity of the base movement
            Eigen::Vector3d vec = v.head(3);
            double dampingForceMag = -100.0 * vec.norm();
            force += dampingForceMag * vec.normalized();
            forceVector.head(3) = force;

            ////Add a force holding the base to upright position
            Eigen::Vector3d goal_vec(0, elementLength, 0.0);
            Eigen::Vector3d error_vec = goal_vec - x.segment<3>(6);
            double error_mag = error_vec.norm();
            error_vec.normalize();
            Eigen::Vector3d righting_force = 50000.0 * error_mag * error_vec;
            forceVector.segment<3>(6) = righting_force;

            //Add a velocity damping force from the swinging motion
            double swing_v_mag = v.segment<3>(6*numElements/2).norm(); //velocity of the midpoint
            Eigen::Vector3d swing_damping =  -100 * swing_v_mag * v.segment<3>(6 * numElements / 2);
            forceVector.segment<3>(6*numElements / 2) += swing_damping;


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
            //forceVector(base + 0) = amp * pointLoad(0);
            //forceVector(base + 1) = amp * endPointLoad(1);
            //forceVector(base + 2) = amp * endPointLoad(2);

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
    beam.simulateTimeDomain(graphics, 300.0, 1/120.0);

    graphics.waitForCompletion();
    graphics.closeGL();

    return 0;
}