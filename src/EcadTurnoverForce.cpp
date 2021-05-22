/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "EcadTurnoverForce.hpp"
#include "EcadTurnoverUtilities.hpp"
template<unsigned DIM>
EcadTurnoverForce<DIM>::EcadTurnoverForce()
   : AbstractForce<DIM>(),
     mAreaElasticityParameter(EcadTurnoverUtilities::Kappa),
     mPerimeterContractilityParameter(EcadTurnoverUtilities::Gamma),
     mLineTensionParameter(EcadTurnoverUtilities::lambda_0),
     mBoundaryLineTensionParameter(EcadTurnoverUtilities::lambda_0)
{
}

template<unsigned DIM>
EcadTurnoverForce<DIM>::~EcadTurnoverForce()
{
}

template<unsigned DIM>
void EcadTurnoverForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    ///\todo: check whether this line influences profiling tests - if so, we should remove it.
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("EcadTurnoverForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    std::vector<double> element_areas_s(num_elements);
    std::vector<double> element_perimeters_s(num_elements);
    std::vector<double> target_areas_s(num_elements);

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);

    for (unsigned int elem_index=0; elem_index<num_elements; ++elem_index)
    {
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
        try
        {
            // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
            // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
            // understandable message. There is a slight chance that the exception is thrown although the error is not about the
            // target areas.
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
            const bool is_L_cell = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("L1 boundary cell")==1;
            if (is_L_cell)
                target_areas[elem_index] *= 1.0;
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use a EcadTurnoverForce");
        }
    }

    std::vector<std::vector<unsigned int> > node_containing_elem_indices(num_nodes);

    for (unsigned int i=0; i<num_nodes; ++i)
    {
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(i)->rGetContainingElementIndices();
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
                iter != containing_elem_indices.end();
                ++iter)
        {
            node_containing_elem_indices[i].push_back(*iter);
        }
    }
    std::vector<c_vector<double, DIM> > area_elasticity_contribution(num_nodes, zero_vector<double>(DIM));
    std::vector<c_vector<double, DIM> > perimeter_contractility_contribution(num_nodes, zero_vector<double>(DIM));
    std::vector<c_vector<double, DIM> > line_tension_contribution(num_nodes, zero_vector<double>(DIM));
#pragma omp parallel for num_threads(12)
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        for (unsigned int elem_number = 0; elem_number<node_containing_elem_indices[node_index].size(); ++elem_number)
        {
            // Get this element, its index and its number of nodes
            const unsigned int elem_index = node_containing_elem_indices[node_index][elem_number];
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(elem_index);
            CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(elem_index);
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's area elasticity (note the minus sign)
            c_vector<double, DIM> element_area_gradient =
                    p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);

            //Boundary cells should not get too small due to contraction
            const double area_elasticity_factor = p_cell->GetCellData()->GetItem("L1 boundary cell")==1.0 ? 1.0 : 1.0;
            area_elasticity_contribution[node_index] -= area_elasticity_factor*GetAreaElasticityParameter()*(element_areas[elem_index] -
                    target_areas[elem_index])*element_area_gradient;

            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;

            // Compute the line tension parameter for each of these edges - be aware that this is half of the actual
            // value for internal edges since we are looping over each of the internal edges twice
            // Here we also take into account boundary line tension parameter
            bool is_prev_edge_on_boundary = p_element->GetEdge(previous_node_local_index)->IsBoundaryEdge();
            bool is_next_edge_on_boundary = p_element->GetEdge(local_index)->IsBoundaryEdge();

            double previous_edge_line_tension_parameter =
                    GetLineTensionParameter(previous_node_local_index, p_cell);
            double next_edge_line_tension_parameter =
                    GetLineTensionParameter(local_index, p_cell);

            // Compute the gradient of each these edges, computed at the present node
            // The first negative  sign due to reversal of the gradient computation
            c_vector<double, DIM> previous_edge_gradient
            =-p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
            c_vector<double, DIM> next_edge_gradient
            = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary line tension (note the minus sign)
            line_tension_contribution[node_index] -= previous_edge_line_tension_parameter*previous_edge_gradient +
                    next_edge_line_tension_parameter*next_edge_gradient;

            // Add the force contribution from this cell's perimeter contractility (note the minus sign)
            c_vector<double, DIM> element_perimeter_gradient;
            element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
            perimeter_contractility_contribution[node_index] -= GetPerimeterContractilityParameter(p_cell)*
                    (element_perimeters[elem_index])*
                    element_perimeter_gradient;
        }
    }
    const c_vector<double, DIM> pop_centroid = rCellPopulation.GetCentroidOfCellPopulation();
    //Forces due to stretching of the tissue
#pragma omp parallel for num_threads(12)
    for (unsigned int node_index=0; node_index<num_nodes; ++node_index)
    {
        c_vector<double, DIM> force_on_node
        = area_elasticity_contribution[node_index] + perimeter_contractility_contribution[node_index] + line_tension_contribution[node_index];
        auto pNode = p_cell_population->GetNode(node_index);
        if (pNode->GetNumNodeAttributes()>0&&pNode->IsBoundaryNode())
        {
            const std::vector<double> attributes = pNode->rGetNodeAttributes();
            const double node_anchorage = attributes[0];
            //Inner nodes. No forces due to stretching added to them
            if (node_anchorage==-2.0)
            {
                force_on_node[0] += 0;
                force_on_node[1] += 0;
            }
            // Boundary nodes
            if (node_anchorage>=0)
            {
                c_vector<double, DIM> equil_position;
                equil_position(0) = attributes[1];
                equil_position(1) = attributes[2];
                c_vector<double, DIM> location= pNode->rGetLocation();
                const c_vector<double, DIM> displacement = location-equil_position;
                if (node_anchorage == 1)
                {
                    // Force that moves right boundary to the right during tissue stretching
                    if (attributes[1]>0.001)
                    {
                        //force_on_node += -EcadTurnoverUtilities::Zone1Constant*displacement;
                        force_on_node[0] = -attributes[1];
                        force_on_node[1] = 0;
                    }
                    else if (attributes[1]<-1.0)
                    {
                        force_on_node[0] = 0;
                        force_on_node[1] = 0;
                    }
                }
                if (node_anchorage == 2)
                {
                    //Do nothing for upper and lower boundary
                }
                if (node_anchorage == 3)
                {
                    // Force that moves right boundary to the right during tissue stretching
                    if (attributes[1]>0.001)
                    {
                        force_on_node[0] = attributes[1];
                        force_on_node[1] = 0;
                    }
                    else if (attributes[1]<-1.0)
                    {
                        force_on_node[0] = 0;
                        force_on_node[1] = 0;
                    }
                }
            }
        }
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
    }
}

template<unsigned DIM>
double EcadTurnoverForce<DIM>::GetAreaElasticityParameter()
{
    return mAreaElasticityParameter;
}

template<unsigned DIM>
double EcadTurnoverForce<DIM>::GetPerimeterContractilityParameter(CellPtr p_cell)
{
    return EcadTurnoverUtilities::Gamma;
}

template<unsigned DIM>
double EcadTurnoverForce<DIM>::GetLineTensionParameter()
{
    return mLineTensionParameter;
}

template<unsigned DIM>
double EcadTurnoverForce<DIM>::GetLineTensionParameter(const unsigned int local_edge_index, CellPtr p_cell)
{
    //Halving is due to summing each edge twice
    const double factor
    = mLineTensionParameter;
    const double grad_tension_param = p_cell->GetCellEdgeData()->GetItemAtIndex("Lambda", local_edge_index);

    return 0.5*grad_tension_param;
}

template<unsigned DIM>
double EcadTurnoverForce<DIM>::GetBoundaryLineTensionParameter()
{
    return mBoundaryLineTensionParameter;
}

template<unsigned DIM>
void EcadTurnoverForce<DIM>::SetAreaElasticityParameter(double areaElasticityParameter)
{
    mAreaElasticityParameter = areaElasticityParameter;
}

template<unsigned DIM>
void EcadTurnoverForce<DIM>::SetPerimeterContractilityParameter(double perimeterContractilityParameter)
{
    mPerimeterContractilityParameter = perimeterContractilityParameter;
}

template<unsigned DIM>
void EcadTurnoverForce<DIM>::SetLineTensionParameter(double lineTensionParameter)
{
    mLineTensionParameter = lineTensionParameter;
}

template<unsigned DIM>
void EcadTurnoverForce<DIM>::SetBoundaryLineTensionParameter(double boundaryLineTensionParameter)
{
    mBoundaryLineTensionParameter = boundaryLineTensionParameter;
}

template<unsigned DIM>
void EcadTurnoverForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<AreaElasticityParameter>" << mAreaElasticityParameter << "</AreaElasticityParameter>\n";
    *rParamsFile << "\t\t\t<PerimeterContractilityParameter>" << mPerimeterContractilityParameter << "</PerimeterContractilityParameter>\n";
    *rParamsFile << "\t\t\t<LineTensionParameter>" << mLineTensionParameter << "</LineTensionParameter>\n";
    *rParamsFile << "\t\t\t<BoundaryLineTensionParameter>" << mBoundaryLineTensionParameter << "</BoundaryLineTensionParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class EcadTurnoverForce<1>;
template class EcadTurnoverForce<2>;
template class EcadTurnoverForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(EcadTurnoverForce)
