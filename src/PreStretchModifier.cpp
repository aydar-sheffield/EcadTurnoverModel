/*
 * PreStretchModifier.cpp
 *
 *  Created on: 21 Nov 2019
 *      Author: aydar
 */

#include "PreStretchModifier.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"

#include "MutableVertexMesh.hpp"
#include "EcadTurnoverUtilities.hpp"
#include <math.h>

template<unsigned DIM>
PreStretchModifier<DIM>::PreStretchModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
PreStretchModifier<DIM>::~PreStretchModifier()
{
}

template<unsigned DIM>
void PreStretchModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Update the cell
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PreStretchModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                                                         std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PreStretchModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    MutableVertexMesh<DIM, DIM>* pMesh
            = static_cast<MutableVertexMesh<DIM, DIM>* >(&(rCellPopulation.rGetMesh()));
    const c_vector<double, DIM> pop_centroid = rCellPopulation.GetCentroidOfCellPopulation();
    const unsigned int n_elements = pMesh->GetNumElements();

    std::vector<double> distances_to_centre(rCellPopulation.GetNumAllCells(),0);
    std::vector<double> cell_boundaries_x;
    unsigned int real_cell_counter = 0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
    {
        const unsigned int location_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        auto element = pMesh->GetElement(location_index);
        const bool is_boundary_cell= element->IsElementOnBoundary();
        if (is_boundary_cell)
        {
            const c_vector<double, DIM> cell_centroid = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            cell_boundaries_x.push_back(cell_centroid(0));
        }
    }

    auto min_element = std::min_element(cell_boundaries_x.begin(), cell_boundaries_x.end());
    auto max_element = std::max_element(cell_boundaries_x.begin(), cell_boundaries_x.end());
    const double min_x = *min_element;
    const double max_x = *max_element;
    const double width = max_x-min_x;

    const double x_displacement = width*0.5;

    for (unsigned int location_index= 0; location_index<n_elements; ++location_index)
    {
        auto pCell = rCellPopulation.GetCellUsingLocationIndex(location_index);
        pCell->GetCellData()->SetItem("target area", 1.0);
        pCell->GetCellData()->SetItem("reference area", 1.0);
        pCell->GetCellData()->SetItem("Displacement", x_displacement);
        pCell->GetCellData()->SetItem("Wing", 1.0);
        const c_vector<double, DIM> cell_centroid = rCellPopulation.GetLocationOfCellCentre(pCell);

        auto element = pMesh->GetElement(location_index);
        bool is_boundary_cell= element->IsElementOnBoundary();
        if (is_boundary_cell)
        {
            const unsigned int n_nodes = element->GetNumNodes();
            c_vector<double, DIM> position_vector = cell_centroid-pop_centroid;
            const double vector_length = norm_2(position_vector);
            position_vector[0] /= vector_length;
            position_vector[1] /= vector_length;
            const double polar_angle = std::acos(position_vector[0]);
            double zone= 0;

            if (polar_angle<EcadTurnoverUtilities::PI*1/4)
                zone = 3.0;
            if (polar_angle>=EcadTurnoverUtilities::PI*1/4&&polar_angle<=EcadTurnoverUtilities::PI*3/4)
                zone = 2.0;
            if (polar_angle>EcadTurnoverUtilities::PI*3/4)
                zone = 1.0;
            pCell->GetCellData()->SetItem("Zone", zone);

            for (unsigned int i=0; i<n_nodes; ++i)
            {
                auto node = element->GetNode(i);
                if (node->IsBoundaryNode())
                {
                    if (node->HasNodeAttributes()&&zone!=1.0)
                        node->rGetNodeAttributes().clear();

                    node->AddNodeAttribute(zone);
                    const c_vector<double, DIM> node_location = node->rGetLocation();
                    node->AddNodeAttribute(node_location(0));
                    node->AddNodeAttribute(node_location(1));
                }
            }
        }
        else
        {
            pCell->GetCellData()->SetItem("Zone", -2.0);
        }

    }

}

template<unsigned DIM>
void PreStretchModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}
// Explicit instantiation
template class PreStretchModifier<1>;
template class PreStretchModifier<2>;
template class PreStretchModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(PreStretchModifier)


