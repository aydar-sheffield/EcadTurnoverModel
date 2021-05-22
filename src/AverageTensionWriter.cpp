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

#include "AverageTensionWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <cmath>

#include "EcadTurnoverUtilities.hpp"
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AverageTensionWriter<ELEMENT_DIM, SPACE_DIM>::AverageTensionWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("Tension.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageTensionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageTensionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageTensionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageTensionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageTensionWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    MutableVertexMesh<SPACE_DIM, SPACE_DIM>* pMesh
            = static_cast<MutableVertexMesh<SPACE_DIM, SPACE_DIM>* >(&(pCellPopulation->rGetMesh()));

    std::vector<double> averages_vector(4,0.0);
    std::vector<unsigned int> n_junctions_vector(4,0);
    std::vector<std::vector<double> > tension_vectors(4, std::vector<double>());
    for (typename VertexBasedCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
            cell_iter != pCellPopulation->End();
            ++cell_iter)
    {
        const unsigned int location_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        auto element = pMesh->GetElement(location_index);

        const bool is_L4_boundary = cell_iter->GetCellData()->GetItem("L4 boundary cell") == 1.0;
        const bool is_mitotic = cell_iter->GetCellData()->GetItem("Mitotic") == 1.0;

        if (!is_L4_boundary&&!is_mitotic)
        {

            std::vector<double> junctional_tension = cell_iter->GetCellEdgeData()->GetItem("Junctional tension");
            const unsigned int n_edges = element->GetNumEdges();
            for (unsigned int i=0; i<n_edges; ++i)
            {
                c_vector<double, SPACE_DIM> AB
                = element->GetEdge(i)->GetNode(1)->rGetLocation() - element->GetEdge(i)->GetNode(0)->rGetLocation();
                const double length=element->GetEdge(i)->rGetLength();
                AB[0] /= length;
                AB[1] /= length;
                if (AB[0]<0)
                    AB[0] *= -1.0;
                const double angle = std::acos(AB[0]);
                if (angle<EcadTurnoverUtilities::PI/6)
                {

                    n_junctions_vector[0]++;
                    averages_vector[0] += junctional_tension[i];
                    tension_vectors[0].push_back(junctional_tension[i]);
                }
                else if (angle>=EcadTurnoverUtilities::PI/6&&angle<=EcadTurnoverUtilities::PI/3)
                {
                    n_junctions_vector[1]++;
                    averages_vector[1] += junctional_tension[i];
                    tension_vectors[1].push_back(junctional_tension[i]);
                }
                else if (angle>EcadTurnoverUtilities::PI/3&&angle<=EcadTurnoverUtilities::PI/2)
                {
                    n_junctions_vector[2]++;
                    averages_vector[2] += junctional_tension[i];
                    tension_vectors[2].push_back(junctional_tension[i]);
                }

                averages_vector[3]+=junctional_tension[i];
                tension_vectors[3].push_back(junctional_tension[i]);
            }

            n_junctions_vector[3] += n_edges;
        }
    }

    std::vector<double> errors(4,0.0);
    for (unsigned int i=0; i<4; ++i)
    {
        if (n_junctions_vector[i]>1)
        {
            averages_vector[i] /= n_junctions_vector[i];
            double sum_squares = 0;
            for (unsigned int j=0; j<n_junctions_vector[i]; ++j)
            {
                sum_squares += std::pow(tension_vectors[i][j] - averages_vector[i],2);
            }
            double std_dev = std::sqrt(sum_squares/(n_junctions_vector[i]-1));
            errors[i] = std_dev/std::sqrt((double)n_junctions_vector[i]);
        }
        else
        {
            averages_vector[i] = 0;
            errors[i] = 0;
            n_junctions_vector[i] = 0;
        }
    }

    for (unsigned int i=0; i<4; ++i)
    {
        *this->mpOutStream << averages_vector[i]<<"\t";
        *this->mpOutStream << errors[i]<<"\t";
        *this->mpOutStream << n_junctions_vector[i]<<"\t";
    }
}

// Explicit instantiation
template class AverageTensionWriter<1,1>;
template class AverageTensionWriter<1,2>;
template class AverageTensionWriter<2,2>;
template class AverageTensionWriter<1,3>;
template class AverageTensionWriter<2,3>;
template class AverageTensionWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AverageTensionWriter)
