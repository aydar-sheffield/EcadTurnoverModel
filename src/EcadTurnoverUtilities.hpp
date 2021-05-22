/*
 * DeltaNotchUtilities.hpp
 *
 *  Created on: 27 Jul 2020
 *      Author: aydar
 */

#ifndef ECADTURNOVERUTILITIES_HPP_
#define ECADTURNOVERUTILITIES_HPP_
#include <vector>
#include <utility>
class EcadTurnoverUtilities
{
public:
    static const double real_to_simulated;
    static const double G1AverageDuration;
    static const double G2Duration;

    static const double friction;
    static const double area_elasticity;
    static const double PI;

    static const double k_minus_0;
    static const double k_plus;
    static const double a;
    static const double k_c_plus;
    static const double k_c_minus;

    static const double mobile_eq;
    static const double bound_eq;

    static const double A0;
    static const double lambda_0;
    static const double lambda;

    static const double Gamma;
    static const double Kappa;

    static const double k_L;
    static const double k_endo;
    static const double boundary_k_L;
    static const double Youngs;
    static const double boundary_Youngs;
    static const double critical_strain;

    static const double Tend;
    static const double k_Lambda;
    static const double dt;

    static const double stretch_time;
    static const double release_time;
    static const double gamma_half_max;
    static const double k_L_mech;

    static double endo_rate(const double gamma);
    static double contractility_param(const double current_length, const double rest_length,
                                      const bool boundary = false);

    static double junctional_tension(const double current_length,const double rest_length,
                                     const bool boundary = false);

    static double rest_length_remodelling_rate(const double current_length, const double rest_length,
                                               const bool is_boundary = false,
                                               const double tension = 0);
    static std::vector<double> ComputePolarEdgeLength(const std::vector<std::vector<double> > cell_nodes,
                                                      const std::vector<double> cell_centroid);
    static std::pair<double, double> ComputeNematicOrderComponents(const std::vector<double> edge_values,
                                                                   const std::vector<double> edge_angles);
    static std::pair<double, double> ComputeAngleAndMagnitude(const std::vector<double> node_0, const std::vector<double> cell_centroid,
                                                              const std::pair<double,double> Q_components);
    static double wrap_angle(const double theta);
private:
    EcadTurnoverUtilities();
};

#endif /* ECADTURNOVERUTILITIES_HPP_ */
