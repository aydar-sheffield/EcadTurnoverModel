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

#include "JunctionalLengthWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
JunctionalLengthWriter<ELEMENT_DIM, SPACE_DIM>::JunctionalLengthWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("JunctionalLengths.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void JunctionalLengthWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void JunctionalLengthWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void JunctionalLengthWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void JunctionalLengthWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void JunctionalLengthWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    MutableVertexMesh<SPACE_DIM, SPACE_DIM>* pMesh
            = static_cast<MutableVertexMesh<SPACE_DIM, SPACE_DIM>* >(&(pCellPopulation->rGetMesh()));
    const unsigned int n_bins = 10;
    const double a = 0, b = 2.5;
    const double h = (b-a)/n_bins;
    std::vector<unsigned int> normalized_lengths(n_bins,0);
    for (typename VertexBasedCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
            cell_iter != pCellPopulation->End();
            ++cell_iter)
    {
        const unsigned int location_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        auto element = pMesh->GetElement(location_index);
        const bool is_L4_boundary = cell_iter->GetCellData()->GetItem("L4 boundary cell") == 1.0;
        const bool is_mitotic = cell_iter->GetCellData()->GetItem("Mitotic") == 1.0;
        const bool is_outlier = cell_iter->GetCellData()->GetItem("Area")<=0.1;
        if (!is_L4_boundary&&!is_mitotic &&!is_outlier)
        {
            const unsigned int n_edges = element->GetNumEdges();
            std::vector<double> junctional_lengths = cell_iter->GetCellEdgeData()->GetItem("Current length");
            const double average_length = cell_iter->GetCellData()->GetItem("Perimeter")/n_edges;
            for (unsigned int i=0; i<n_edges; ++i)
            {
                junctional_lengths[i] /=average_length;
                unsigned int bin = std::floor(junctional_lengths[i]/h);
                if (junctional_lengths[i]>=b)
                {
                    normalized_lengths[9]++;
                }
                else
                {
                    normalized_lengths[bin]++;
                }
            }

        }
    }
    for (unsigned int i=0; i<n_bins; ++i)
    {
        *this->mpOutStream << normalized_lengths[i]<<"\t";
    }
}

// Explicit instantiation
template class JunctionalLengthWriter<1,1>;
template class JunctionalLengthWriter<1,2>;
template class JunctionalLengthWriter<2,2>;
template class JunctionalLengthWriter<1,3>;
template class JunctionalLengthWriter<2,3>;
template class JunctionalLengthWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(JunctionalLengthWriter)
