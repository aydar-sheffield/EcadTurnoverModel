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

#include "AverageElongationWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AverageElongationWriter<ELEMENT_DIM, SPACE_DIM>::AverageElongationWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("Elongation.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageElongationWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageElongationWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageElongationWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageElongationWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AverageElongationWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    MutableVertexMesh<SPACE_DIM, SPACE_DIM>* pMesh
            = static_cast<MutableVertexMesh<SPACE_DIM, SPACE_DIM>* >(&(pCellPopulation->rGetMesh()));

    double total_elongation=0, average_elongation=0;
    unsigned int n_non_boundary_cells = 0;
    std::vector<double> elongation_vector;
    double error = 0;

    double total_shape_param=0;
    std::vector<double> shape_param_vector;
    double error_shape_param = 0;

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
            const double elongation = cell_iter->GetCellData()->GetItem("Elongation");
            n_non_boundary_cells++;
            total_elongation += elongation;
            elongation_vector.push_back(elongation);

            const double shape_param_value = cell_iter->GetCellData()->GetItem("Perimeter")/std::sqrt(cell_iter->GetCellData()->GetItem("Area"));
            total_shape_param+=shape_param_value;
            shape_param_vector.push_back(shape_param_value);

        }
    }

    average_elongation = total_elongation/n_non_boundary_cells;
    double average_shape_param = total_shape_param/n_non_boundary_cells;
    double sum_squares = 0;
    double std_dev = 0;

    for (unsigned int i=0; i<n_non_boundary_cells; ++i)
    {
        sum_squares += std::pow(elongation_vector[i]-average_elongation,2);
    }
    std_dev = std::sqrt(sum_squares/(n_non_boundary_cells-1));
    error = std_dev/std::sqrt((double)n_non_boundary_cells);

    sum_squares = 0;
    for (unsigned int i=0; i<n_non_boundary_cells; ++i)
    {
        sum_squares += std::pow(shape_param_vector[i]-average_shape_param,2);
    }
    std_dev = std::sqrt(sum_squares/(n_non_boundary_cells-1));
    error_shape_param = std_dev/std::sqrt((double)n_non_boundary_cells);


    *this->mpOutStream << average_elongation<<"\t";
    *this->mpOutStream << error<<"\t";
    *this->mpOutStream << n_non_boundary_cells<<"\t";

    *this->mpOutStream << average_shape_param<<"\t";
    *this->mpOutStream << error_shape_param<<"\t";
    *this->mpOutStream << n_non_boundary_cells<<"\t";
}

// Explicit instantiation
template class AverageElongationWriter<1,1>;
template class AverageElongationWriter<1,2>;
template class AverageElongationWriter<2,2>;
template class AverageElongationWriter<1,3>;
template class AverageElongationWriter<2,3>;
template class AverageElongationWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(AverageElongationWriter)
