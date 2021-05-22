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

#include "AreaAndBoundaryWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AreaAndBoundaryWriter<ELEMENT_DIM, SPACE_DIM>::AreaAndBoundaryWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("AreaAndBoundary.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AreaAndBoundaryWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AreaAndBoundaryWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AreaAndBoundaryWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AreaAndBoundaryWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AreaAndBoundaryWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    MutableVertexMesh<SPACE_DIM, SPACE_DIM>* pMesh
            = static_cast<MutableVertexMesh<SPACE_DIM, SPACE_DIM>* >(&(pCellPopulation->rGetMesh()));
    double total_area = 0;

    std::vector<double> boundary_x;
    std::vector<double> boundary_y;
    for (typename VertexBasedCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
            cell_iter != pCellPopulation->End();
            ++cell_iter)
    {
        const unsigned int location_index = pCellPopulation->GetLocationIndexUsingCell(*cell_iter);
        auto element = pMesh->GetElement(location_index);
        const double area = pCellPopulation->GetVolumeOfCell(*cell_iter);
        total_area += area;
        const bool is_boundary_cell = element->IsElementOnBoundary();
        if (is_boundary_cell)
        {
            c_vector<double, SPACE_DIM> centroid = pCellPopulation->GetLocationOfCellCentre(*cell_iter);
            boundary_x.push_back(centroid(0));
            boundary_y.push_back(centroid(1));
        }

    }
    *this->mpOutStream << total_area <<"\t";
    *this->mpOutStream << boundary_x.size()<<"\t";
    for (unsigned int i=0; i<boundary_x.size(); ++i)
    {
        *this->mpOutStream << boundary_x[i]<<"\t";
        *this->mpOutStream << boundary_y[i]<<"\t";
    }

}

// Explicit instantiation
template class AreaAndBoundaryWriter<1,1>;
template class AreaAndBoundaryWriter<1,2>;
template class AreaAndBoundaryWriter<2,2>;
template class AreaAndBoundaryWriter<1,3>;
template class AreaAndBoundaryWriter<2,3>;
template class AreaAndBoundaryWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AreaAndBoundaryWriter)
