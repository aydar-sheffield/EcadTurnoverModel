/*
 *
 *
 *  Created on: 27 Jul 2020
 *      Author: aydar
 */

#include "EcadTurnoverUtilities.hpp"
#include <cmath>

const double EcadTurnoverUtilities::real_to_simulated = 50;// = 1hr in real time
const double EcadTurnoverUtilities::G1AverageDuration = 3*real_to_simulated;//4-5 hr
const double EcadTurnoverUtilities::G2Duration = EcadTurnoverUtilities::G1AverageDuration*0.5;

//Friction sets the time scale of vertex models 1
const double EcadTurnoverUtilities::friction = 1;
const double EcadTurnoverUtilities::PI = 3.141592653589793;

const double EcadTurnoverUtilities::k_minus_0 = 0.0027;//3.685e-5s in real time ==
//const double EcadTurnoverUtilities::k_minus_0 = 1.0/(8.0*real_to_simulated);//~8 hours

const double EcadTurnoverUtilities::k_plus = 4.525*k_minus_0;

const double EcadTurnoverUtilities::a = 2.1381*k_minus_0;
const double EcadTurnoverUtilities::k_c_minus = 1.0/(215/3600.0*real_to_simulated)-k_plus;//~t_minus ~= 215s
const double EcadTurnoverUtilities::k_c_plus = k_c_minus*2;

const double EcadTurnoverUtilities::mobile_eq = k_c_plus/(k_c_minus+k_plus);
const double EcadTurnoverUtilities::bound_eq = k_plus/k_minus_0*k_c_plus/(k_c_minus+k_plus);

const double EcadTurnoverUtilities::A0 = 1;
const double EcadTurnoverUtilities::lambda_0 = 0;
const double EcadTurnoverUtilities::lambda = 0.01;//max value is lambda_0/(2*bound_eq) ~= 0.0075

const double EcadTurnoverUtilities::Gamma = 0.04;
const double EcadTurnoverUtilities::Kappa = 1;

const double EcadTurnoverUtilities::k_L = 1.0/(8*real_to_simulated);
const double EcadTurnoverUtilities::k_endo = 1.0/(8*real_to_simulated);
const double EcadTurnoverUtilities::boundary_k_L = k_L;
const double EcadTurnoverUtilities::Youngs = 30;//5, 15 -- original, 25, 35
const double EcadTurnoverUtilities::boundary_Youngs = 20;
const double EcadTurnoverUtilities::critical_strain =  0.0;

const double EcadTurnoverUtilities::Tend = 16*real_to_simulated;
const double EcadTurnoverUtilities::k_Lambda = 1.0/(1*real_to_simulated);
const double EcadTurnoverUtilities::dt = 0.01;

const double EcadTurnoverUtilities::stretch_time = real_to_simulated/10;
const double EcadTurnoverUtilities::release_time = real_to_simulated*3.5;
const double EcadTurnoverUtilities::gamma_half_max = 0.75;
const double EcadTurnoverUtilities::k_L_mech = 32*k_L;//32 before

double EcadTurnoverUtilities::endo_rate(const double gamma)
{
    return k_L+gamma/(gamma+gamma_half_max)*k_L_mech;
}


double EcadTurnoverUtilities::contractility_param(const double current_length, const double rest_length,
                                                  const bool boundary)
{
    const double modulus = boundary? boundary_Youngs : Youngs;
    return current_length>rest_length ? modulus/rest_length*std::pow(current_length-rest_length,2) : 0;
}

double EcadTurnoverUtilities::junctional_tension(const double current_length, const double rest_length,
                                                 const bool boundary)
{
    /*const double strain = current_length-rest_length;
    const double modulus = boundary? boundary_Youngs : Youngs;
    double contractility_contribution = strain>0 ? strain*modulus : 0;
    const double line_tension = boundary? std::abs(lambda_0):lambda_0;
    return line_tension - adhesion + contractility_contribution;*/
    const double displacement = current_length-rest_length;

    const double modulus = boundary? boundary_Youngs : Youngs;

    double line_tension =0;

    if (boundary)
    {
        const double elastic = displacement>0 ? modulus*displacement : 0;
        line_tension = std::abs(lambda_0) + elastic;
    }
    else
    {
        const double strain = displacement/rest_length;
        const double contractility = strain>0&&current_length>0.15 ? modulus*(current_length-rest_length) : 0;
        line_tension = lambda_0 + contractility;
    }
    const double strain = displacement/rest_length;
    const double contractility = strain>0&&current_length>0.15 ? modulus*(current_length-rest_length) : 0;
    line_tension = lambda_0 + contractility;

    return line_tension;
}


double EcadTurnoverUtilities::rest_length_remodelling_rate(const double current_length, const double rest_length,
                                                           const bool is_boundary,const double tension)
{
    /*if (is_boundary)
    {
        //return current_length<rest_length? k_L*0 : 0.0;
        const double result = k_L+tension/(tension+gamma_half_max)*k_L_mech;
        return result;
    }
    else
    {
        const double result = k_L+tension/(tension+gamma_half_max)*k_L_mech;
        return result;
    }*/
    //const double result = k_L+tension/(tension+gamma_half_max)*k_L_mech;
    //Test 1/10min, 1/15min, 1/30 min, 1/60min
    const double result = 1.0/(real_to_simulated/4);
    return result;
}

std::vector<double> EcadTurnoverUtilities::ComputePolarEdgeLength(const std::vector<std::vector<double> > node_locations,
                                                                  const std::vector<double> cell_centroid)
{
    const unsigned int n_nodes = node_locations.size();
    std::vector<double> result(n_nodes);
    double OAx, OAy, OBx, OBy;
    for (unsigned int i=0; i<n_nodes-1; ++i)
    {
        OAx = node_locations[i][0] - cell_centroid[0];
        OAy = node_locations[i][1] - cell_centroid[1];
        OBx = node_locations[i+1][0] - cell_centroid[0];
        OBy = node_locations[i+1][1] - cell_centroid[1];
        const double product = OAx*OBx+OAy*OBy;
        const double normOA = std::sqrt(std::pow(OAx,2)+std::pow(OAy,2));
        const double normOB = std::sqrt(std::pow(OBx,2)+std::pow(OBy,2));
        result[i] = std::acos(product/(normOA*normOB));
    }
    OAx = node_locations[n_nodes-1][0] - cell_centroid[0];
    OAy = node_locations[n_nodes-1][1] - cell_centroid[1];
    OBx = node_locations[0][0] - cell_centroid[0];
    OBy = node_locations[0][1] - cell_centroid[1];
    const double product = OAx*OBx+OAy*OBy;
    const double normOA = std::sqrt(std::pow(OAx,2)+std::pow(OAy,2));
    const double normOB = std::sqrt(std::pow(OBx,2)+std::pow(OBy,2));
    result[n_nodes-1] = std::acos(product/(normOA*normOB));
    return result;
}

std::pair<double, double> EcadTurnoverUtilities::ComputeNematicOrderComponents(const std::vector<double> edge_values,
                                                                         const std::vector<double> edge_angles)
{
    std::pair<double, double> result;
    const unsigned int n_edges = edge_values.size();
    std::vector<std::pair<double,double> > edge_polar_location(n_edges);
    edge_polar_location[0].first = 0;
    edge_polar_location[0].second = edge_angles[0];
    for (unsigned int i=1; i<n_edges; ++i)
    {
        edge_polar_location[i].first = edge_polar_location[i-1].second;
        edge_polar_location[i].second = edge_polar_location[i].first+edge_angles[i];
    }

    double Q_1, Q_2;
    for (unsigned int i=0; i<n_edges; ++i)
    {
        Q_1 += edge_values[i]*0.5*(std::sin(2*edge_polar_location[i].second)
        -std::sin(2*edge_polar_location[i].first));
        Q_2 += -edge_values[i]*0.5*(std::cos(2*edge_polar_location[i].second)
        -std::cos(2*edge_polar_location[i].first));
    }
    result.first = Q_1;
    result.second = Q_2;
    return result;
}

std::pair<double, double> EcadTurnoverUtilities::ComputeAngleAndMagnitude(const std::vector<double> node_0, const std::vector<double> cell_centroid,
                                                                          const std::pair<double,double> Q_components)
{
    std::pair<double, double> result;
    const double Q_magnitude = std::sqrt(std::pow(Q_components.first,2)
    +std::pow(Q_components.second,2));
    double polarity_angle = wrap_angle(std::atan2(Q_components.second, Q_components.first))/2;
    const double OAx = node_0[0]-cell_centroid[0];
    const double OAy = node_0[1]-cell_centroid[1];
    const double norm_OA = std::sqrt(std::pow(OAx,2)+std::pow(OAy,2));
    const double node_0_angle = std::acos(OAx/norm_OA);
    polarity_angle = wrap_angle(OAy>0 ? polarity_angle+node_0_angle : 2*PI-polarity_angle+node_0_angle);
    result.first = polarity_angle;
    result.second = Q_magnitude;
    return result;
}

double EcadTurnoverUtilities::wrap_angle(const double theta)
{
    if (theta<0)
        return 2*PI+theta;
    else if(theta>2*PI)
        return theta-2*PI;
    else
        return theta;

}

