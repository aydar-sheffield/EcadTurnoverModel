/*
 * EcadTurnoverModifier.cpp
 *
 *  Created on: 21 Nov 2019
 *      Author: aydar
 */

#include "EcadTurnoverModifier.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"

#include "MutableVertexMesh.hpp"
#include "CellSrnModel.hpp"
#include "EcadJunctionSrnModel.hpp"
#include "EcadTurnoverUtilities.hpp"
#include <math.h>
#include <omp.h>
#include <cstdint>
template<unsigned DIM>
EcadTurnoverModifier<DIM>::EcadTurnoverModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
EcadTurnoverModifier<DIM>::~EcadTurnoverModifier()
{
}

template<unsigned DIM>
void EcadTurnoverModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Update the cell
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void EcadTurnoverModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                         std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void EcadTurnoverModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

    const double current_time = SimulationTime::Instance()->GetTime();
    std::cout<<"Time: "<<current_time<<std::endl;
    MutableVertexMesh<DIM, DIM>* pMesh
            = static_cast<MutableVertexMesh<DIM, DIM>* >(&(rCellPopulation.rGetMesh()));
    const c_vector<double, DIM> pop_centroid = rCellPopulation.GetCentroidOfCellPopulation();
    const unsigned int n_elements = pMesh->GetNumElements();

    std::vector<double> distances_to_centre(rCellPopulation.GetNumAllCells(),0);
    std::vector<double> cell_boundaries_x;

    for (unsigned int location_index= 0; location_index<n_elements; ++location_index)
    {
        auto cell_iter = rCellPopulation.GetCellUsingLocationIndex(location_index);
        auto element = pMesh->GetElement(location_index);
        const bool is_boundary_cell= element->IsElementOnBoundary();

        const c_vector<double, DIM> cell_centroid = rCellPopulation.GetLocationOfCellCentre(cell_iter);
        if (is_boundary_cell)
        {
            cell_boundaries_x.push_back(cell_centroid(0));
        }

        auto p_cell_edge_model = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());

        /* Cell junction Ecads */
        std::vector<double> mobile_ecad;
        std::vector<double> bound_ecad;
        /* Cell junction Rest lengths*/
        std::vector<double> edge_lengths(element->GetNumEdges(), 0.0);
        std::vector<double> rest_lengths = edge_lengths;
        std::vector<double> mobile_fraction(rest_lengths);

        for (unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            boost::shared_ptr<EcadJunctionSrnModel> p_model
            = boost::static_pointer_cast<EcadJunctionSrnModel>(p_cell_edge_model->GetEdgeSrn(i));
            double mobile = p_model->GetMobileEcad();
            double bound = p_model->GetBoundEcad();

            mobile_ecad.push_back(mobile);
            bound_ecad.push_back(bound);
            mobile_fraction.push_back(mobile/(mobile+bound));

            edge_lengths[i] = element->GetEdge(i)->rGetLength();
            rest_lengths[i] = p_model->GetRestLength();
            if (rest_lengths[i]<0.0)
            {
                rest_lengths[i] = 1.0*edge_lengths[i];
                p_model->SetRestLength(rest_lengths[i]);
            }
        }
        cell_iter->GetCellEdgeData()->SetItem("Mobile Ecad", mobile_ecad);
        cell_iter->GetCellEdgeData()->SetItem("Bound Ecad", bound_ecad);
        cell_iter->GetCellEdgeData()->SetItem("Mobile fraction", mobile_fraction);

        cell_iter->GetCellEdgeData()->SetItem("Rest length", rest_lengths);
        cell_iter->GetCellEdgeData()->SetItem("Current length", edge_lengths);

        const double perimeter = pMesh->GetSurfaceAreaOfElement(location_index);
        cell_iter->GetCellData()->SetItem("Perimeter", perimeter);

        //If a new boundary node appears on cells being pulled
        if (is_boundary_cell)
        {
            const unsigned int n_nodes = element->GetNumNodes();
            const double zone = cell_iter->GetCellData()->GetItem("Zone");
            assert(zone>0);
            for (unsigned int i=0; i<n_nodes; ++i)
            {
                auto node = element->GetNode(i);
                if (node->IsBoundaryNode()&&!node->HasNodeAttributes())
                {
                    Node<DIM>* attribute_node = nullptr;
                    Node<DIM>* next_node = element->GetNode((i+1)%n_nodes);
                    Node<DIM>* previous_node = element->GetNode((i-1+n_nodes)%n_nodes);
                    if (next_node->IsBoundaryNode()&&next_node->HasNodeAttributes())
                        attribute_node = next_node;
                    else if (previous_node->IsBoundaryNode()&&previous_node->HasNodeAttributes())
                        attribute_node = previous_node;
                    else
                        EXCEPTION("A new node is created with neighbouring non-boundary node. This shouldn't occur");
                    assert(attribute_node);
                    node->AddNodeAttribute(attribute_node->rGetNodeAttributes()[0]);
                    node->AddNodeAttribute(attribute_node->rGetNodeAttributes()[1]);
                    node->AddNodeAttribute(attribute_node->rGetNodeAttributes()[2]);
                }
            }
            const double stretch_time = EcadTurnoverUtilities::stretch_time;
            const double release_time = EcadTurnoverUtilities::release_time;
            if (zone == 1.0)
            {
                for (unsigned int i=0; i<n_nodes; ++i)
                {
                    auto pNode = element->GetNode(i);
                    if (pNode->GetNumNodeAttributes()>0&&pNode->IsBoundaryNode())
                    {
                        std::vector<double> attributes = pNode->rGetNodeAttributes();
                        const double node_anchorage = attributes[0];

                        if (node_anchorage==1)
                        {
                            const double d = cell_iter->GetCellData()->GetItem("Displacement");
                            /*double velocity
                            = d*2/EcadTurnoverUtilities::Tend*(1.0-current_time/EcadTurnoverUtilities::Tend);*/
                            double velocity = d/stretch_time;
                            if (current_time>=stretch_time&&current_time<=release_time)
                            {
                                velocity =-2.0;
                            }
                            else if (current_time>=release_time)
                            {
                                velocity = 0.0;
                            }
                            pNode->rGetNodeAttributes()[1] = EcadTurnoverUtilities::friction*velocity;
                        }
                    }
                }
            }

            if (zone == 3.0)
            {
                for (unsigned int i=0; i<n_nodes; ++i)
                {
                    auto pNode = element->GetNode(i);
                    if (pNode->GetNumNodeAttributes()>0&&pNode->IsBoundaryNode())
                    {
                        std::vector<double> attributes = pNode->rGetNodeAttributes();
                        const double node_anchorage = attributes[0];

                        if (node_anchorage==3)
                        {
                            const double d = cell_iter->GetCellData()->GetItem("Displacement");
                            /*const double velocity
                            = d*2/EcadTurnoverUtilities::Tend*(1.0-current_time/EcadTurnoverUtilities::Tend);*/
                            double velocity = d/(stretch_time);
                            if (current_time>=stretch_time&&current_time<=release_time)
                            {
                                velocity =-2.0;
                            }
                            else if (current_time>=release_time)
                            {
                                velocity = 0.0;
                            }
                            pNode->rGetNodeAttributes()[1] = EcadTurnoverUtilities::friction*velocity;
                        }
                    }
                }
            }

        }

    }


#pragma omp parallel for num_threads(12)
    for (unsigned int location_index= 0; location_index<n_elements; ++location_index)
    {
        auto cell_iter = rCellPopulation.GetCellUsingLocationIndex(location_index);
        auto element = pMesh->GetElement(location_index);
        auto p_cell_edge_model = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());
        const unsigned int n_cell_edges = p_cell_edge_model->GetNumEdgeSrn();

        std::vector<double> current_lengths = cell_iter->GetCellEdgeData()->GetItem("Current length");
        std::vector<double> rest_lengths = cell_iter->GetCellEdgeData()->GetItem("Rest length");
        std::vector<double> junction_ecad = cell_iter->GetCellEdgeData()->GetItem("Bound Ecad");
        //const double cell_perim_deviation = cell_iter->GetCellData()->GetItem("Perimeter deviation");
        const double perimeter = cell_iter->GetCellData()->GetItem("Perimeter");
        std::vector<double> neighbour_ecad(n_cell_edges);
        std::vector<double> neighbour_perim_deviations(n_cell_edges);

        std::vector<double> junctional_tension(n_cell_edges);
        std::vector<double> endocytosis_rate(n_cell_edges);
        std::vector<double> junctional_tension_param(n_cell_edges);
        std::vector<double> perimeter_contractility(n_cell_edges);
        std::vector<double> remodelling_rate(n_cell_edges);
        bool is_neighbour_boundary = false;
        bool is_boundary_cell = element->IsElementOnBoundary();

        const double elongation
        = is_boundary_cell|| is_neighbour_boundary? 0 : pMesh->GetElongationShapeFactorOfElement(location_index);
        const c_vector<double, DIM> short_axis = pMesh->GetShortAxisOfElement(location_index);
        c_vector<double, DIM> long_axis;
        long_axis(0) = -short_axis(1);
        long_axis(1) = short_axis(0);
        if (long_axis(1)<0)
        {
            long_axis(0) *=-1.0;
            long_axis(1) *=-1.0;
        }
        cell_iter->GetCellData()->SetItem("Long axis x", long_axis(0));
        cell_iter->GetCellData()->SetItem("Long axis y", long_axis(1));

        cell_iter->GetCellData()->SetItem("Elongation", elongation);
        for (unsigned int i=0; i<n_cell_edges; ++i)
        {
            bool is_edge_neighbour_cell_on_boundary = false;
            //Get neighbouring cell's values of delta on this
            auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(cell_iter, i);
            double neigh_ecad_value = 0;
            double neigh_per_value = 0;
            for (auto neighbourIndex: elemNeighbours)
            {
                auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);
                std::vector<double> neighbour_bound_ecad = neighbourCell->GetCellEdgeData()->GetItem("Bound Ecad");
                neigh_ecad_value += neighbour_bound_ecad[neighbourIndex.second];

                neigh_per_value = neighbourCell->GetCellData()->GetItem("Perimeter");
                if (pMesh->GetElement(neighbourIndex.first)->IsElementOnBoundary())
                {
                    is_neighbour_boundary = true;
                    is_edge_neighbour_cell_on_boundary = true;
                }
            }
            bool is_boundary_edge = element->GetEdge(i)->IsBoundaryEdge();
            if (is_boundary_edge)
            {
                neigh_ecad_value = junction_ecad[i];
            }
            neighbour_ecad[i] = neigh_ecad_value;
            assert(rest_lengths[i]>0);
            const bool is_boundary = is_boundary_cell||is_edge_neighbour_cell_on_boundary;
            junctional_tension[i] = EcadTurnoverUtilities::junctional_tension(current_lengths[i], rest_lengths[i],
                                                                              is_boundary);
            junctional_tension_param[i] = junctional_tension[i];
            perimeter_contractility[i]
                                    = EcadTurnoverUtilities::Gamma*(perimeter);
            junctional_tension[i]=0.5*junctional_tension_param[i]+perimeter_contractility[i];
            endocytosis_rate[i]= is_boundary ? EcadTurnoverUtilities::endo_rate(1.0) :
                    EcadTurnoverUtilities::endo_rate(junctional_tension[i]);
            const double remodel_rate = EcadTurnoverUtilities::rest_length_remodelling_rate(current_lengths[i], rest_lengths[i],
                                                                                            is_boundary,junctional_tension[i]);
            remodelling_rate[i] = remodel_rate;
        }
        cell_iter->GetCellEdgeData()->SetItem("Remodelling rate", remodelling_rate);

        const double area = pMesh->GetVolumeOfElement(location_index);
        cell_iter->GetCellData()->SetItem("Area", area);
        cell_iter->GetCellData()->SetItem("L1 boundary cell", is_boundary_cell);
        cell_iter->GetCellData()->SetItem("L2 boundary cell", is_neighbour_boundary);
        const bool is_mitotic = cell_iter->GetCellData()->GetItem("target area")>EcadTurnoverUtilities::A0;
        cell_iter->GetCellData()->SetItem("Mitotic", is_mitotic);
        cell_iter->GetCellEdgeData()->SetItem("L1 boundary edge", std::vector<double>(n_cell_edges, is_boundary_cell));
        cell_iter->GetCellEdgeData()->SetItem("L2 boundary edge", std::vector<double>(n_cell_edges, is_neighbour_boundary));

        cell_iter->GetCellEdgeData()->SetItem("Junctional tension", junctional_tension);
        cell_iter->GetCellEdgeData()->SetItem("Endocytosis rate", endocytosis_rate);

        cell_iter->GetCellEdgeData()->SetItem("Lambda", junctional_tension_param);
        cell_iter->GetCellEdgeData()->SetItem("Perimeter contractility", perimeter_contractility);

        double average_lambda = 0;
        double average_tension = 0;
        for (unsigned int i=0; i<n_cell_edges; ++i)
        {
            average_lambda += junctional_tension_param[i];
            average_tension += junctional_tension[i];
        }
        average_lambda /= n_cell_edges;
        average_tension /= n_cell_edges;

        cell_iter->GetCellData()->SetItem("Average tension", average_tension);

        //Polarity
        double Q_magnitude= 0, nematic_x= 0, nematic_y= 0, polarity_angle= 0;
        //Transform to polar coordinates...
        auto centroid_c = rCellPopulation.GetLocationOfCellCentre(cell_iter);
        const std::vector<double> c_centroid = {centroid_c[0], centroid_c[1]};
        std::vector<std::vector<double> > node_positions(element->GetNumNodes(),
                                                         std::vector<double>(2,0.0));
        for (unsigned int i=0; i<node_positions.size(); ++i)
        {
            auto cell_node = element->GetNode(i);
            node_positions[i][0] = cell_node->rGetLocation()[0];
            node_positions[i][1] = cell_node->rGetLocation()[1];
        }
        std::vector<double> edge_angles = EcadTurnoverUtilities::ComputePolarEdgeLength(node_positions, c_centroid);

        std::vector<double> bound_ecad = cell_iter->GetCellEdgeData()->GetItem("Bound Ecad");
        std::vector<double> mobile_ecad = cell_iter->GetCellEdgeData()->GetItem("Mobile Ecad");
        std::vector<double> total_ecad(mobile_ecad);
        for (unsigned int i=0; i<total_ecad.size(); ++i)
        {
            total_ecad[i] = bound_ecad[i]+mobile_ecad[i];
        }
        if (!is_boundary_cell)
        {
            std::pair<double, double> angle_magnitude(0.0,0.0), angle_magnitude_raw(0.0,0.0);
            std::pair<double, double> Q_components = EcadTurnoverUtilities::ComputeNematicOrderComponents(total_ecad,
                                                                                                          edge_angles);
            angle_magnitude = EcadTurnoverUtilities::ComputeAngleAndMagnitude(node_positions[0], c_centroid, Q_components);

            Q_magnitude = angle_magnitude.second;
            polarity_angle = angle_magnitude.first;
            if (polarity_angle<0)
                polarity_angle += 2.0*EcadTurnoverUtilities::PI;
            nematic_x = std::cos(polarity_angle)*Q_magnitude;
            nematic_y = std::sin(polarity_angle)*Q_magnitude;

        }
        cell_iter->GetCellData()->SetItem("Polarity magnitude", Q_magnitude);
        cell_iter->GetCellData()->SetItem("Nematic x", nematic_x);
        cell_iter->GetCellData()->SetItem("Nematic y", nematic_y);
        cell_iter->GetCellData()->SetItem("Polarity angle", polarity_angle);
    }
#pragma omp parallel for num_threads(12)
    for (unsigned int location_index= 0; location_index<n_elements; ++location_index)
    {
        auto cell_iter = rCellPopulation.GetCellUsingLocationIndex(location_index);
        auto element = pMesh->GetElement(location_index);
        const unsigned int n_cell_edges = element->GetNumEdges();
        bool is_L3_cell = false;
        if (cell_iter->GetCellData()->GetItem("L2 boundary cell")==0.0)
        {
            for (unsigned int i=0; i<n_cell_edges; ++i)
            {
                //Get neighbouring cell's values of delta on this
                auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(cell_iter, i);
                for (auto neighbourIndex: elemNeighbours)
                {
                    auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);
                    if (neighbourCell->GetCellData()->GetItem("L2 boundary cell") == 1.0)
                        is_L3_cell = true;
                }
            }
        }
        else
        {
            is_L3_cell = true;
        }
        cell_iter->GetCellData()->SetItem("L3 boundary cell", is_L3_cell);
    }

#pragma omp parallel for num_threads(12)
    for (unsigned int location_index= 0; location_index<n_elements; ++location_index)
    {
        auto cell_iter = rCellPopulation.GetCellUsingLocationIndex(location_index);
        auto element = pMesh->GetElement(location_index);
        const unsigned int n_cell_edges = element->GetNumEdges();
        bool is_L4_cell = false;
        if (cell_iter->GetCellData()->GetItem("L3 boundary cell")==0.0)
        {
            for (unsigned int i=0; i<n_cell_edges; ++i)
            {
                //Get neighbouring cell's values of delta on this
                auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(cell_iter, i);
                for (auto neighbourIndex: elemNeighbours)
                {
                    auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);
                    if (neighbourCell->GetCellData()->GetItem("L3 boundary cell") == 1.0)
                        is_L4_cell = true;
                }
            }
        }
        else
        {
            is_L4_cell = true;
        }
        cell_iter->GetCellData()->SetItem("L4 boundary cell", is_L4_cell);
        /*double ref_area = current_time<=shrink_duration ? 0.5-0.36*current_time/(EcadTurnoverUtilities::Tend/4):0.27;
        cell_iter->GetCellData()->SetItem("target area", is_L4_cell ? 1.0: ref_area);*/
        const std::vector<double> is_L4_edge(n_cell_edges, is_L4_cell);
        cell_iter->GetCellEdgeData()->SetItem("L4 edge", is_L4_edge);
    }

#pragma omp parallel for num_threads(12)
    for (unsigned int location_index= 0; location_index<n_elements; ++location_index)
    {
        auto cell_iter = rCellPopulation.GetCellUsingLocationIndex(location_index);
        auto element = pMesh->GetElement(location_index);
        const unsigned int n_cell_edges = element->GetNumEdges();
        bool is_L5_cell = false;
        if (cell_iter->GetCellData()->GetItem("L4 boundary cell")==0.0)
        {
            for (unsigned int i=0; i<n_cell_edges; ++i)
            {
                //Get neighbouring cell's values of delta on this
                auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(cell_iter, i);
                for (auto neighbourIndex: elemNeighbours)
                {
                    auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);
                    if (neighbourCell->GetCellData()->GetItem("L4 boundary cell") == 1.0)
                        is_L5_cell = true;
                }
            }
        }
        else
        {
            is_L5_cell = true;
        }
        cell_iter->GetCellData()->SetItem("L5 boundary cell", is_L5_cell);
        /*double ref_area = current_time<=shrink_duration ? 0.5-0.36*current_time/(EcadTurnoverUtilities::Tend/4):0.27;
        cell_iter->GetCellData()->SetItem("target area", is_L4_cell ? 1.0: ref_area);*/
        const std::vector<double> is_L4_edge(n_cell_edges, is_L5_cell);
        cell_iter->GetCellEdgeData()->SetItem("L5 edge", is_L4_edge);
    }
#pragma omp parallel for num_threads(12)
    for (unsigned int location_index= 0; location_index<n_elements; ++location_index)
    {
        auto cell_iter = rCellPopulation.GetCellUsingLocationIndex(location_index);
        auto element = pMesh->GetElement(location_index);
        const unsigned int n_cell_edges = element->GetNumEdges();
        bool is_L6_cell = false;
        if (cell_iter->GetCellData()->GetItem("L5 boundary cell")==0.0)
        {
            for (unsigned int i=0; i<n_cell_edges; ++i)
            {
                //Get neighbouring cell's values of delta on this
                auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(cell_iter, i);
                for (auto neighbourIndex: elemNeighbours)
                {
                    auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);
                    if (neighbourCell->GetCellData()->GetItem("L5 boundary cell") == 1.0)
                        is_L6_cell = true;
                }
            }
        }
        else
        {
            is_L6_cell = true;
        }
        cell_iter->GetCellData()->SetItem("L6 boundary cell", is_L6_cell);
        /*double ref_area = current_time<=shrink_duration ? 0.5-0.36*current_time/(EcadTurnoverUtilities::Tend/4):0.27;
        cell_iter->GetCellData()->SetItem("target area", is_L4_cell ? 1.0: ref_area);*/
        const std::vector<double> is_L6_edge(n_cell_edges, is_L6_cell);
        cell_iter->GetCellEdgeData()->SetItem("L6 edge", is_L6_edge);
    }
}

template<unsigned DIM>
void EcadTurnoverModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}
// Explicit instantiation
template class EcadTurnoverModifier<1>;
template class EcadTurnoverModifier<2>;
template class EcadTurnoverModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(EcadTurnoverModifier)


