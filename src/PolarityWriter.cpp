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

#include "PolarityWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "EcadTurnoverUtilities.hpp"
#include <cmath>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PolarityWriter<ELEMENT_DIM, SPACE_DIM>::PolarityWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("Polarity.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PolarityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PolarityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PolarityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PolarityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void PolarityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    MutableVertexMesh<SPACE_DIM, SPACE_DIM>* pMesh
            = static_cast<MutableVertexMesh<SPACE_DIM, SPACE_DIM>* >(&(pCellPopulation->rGetMesh()));
    double mean_magnitude = 0;
    double mean_angle = 0;
    unsigned int num_cells = 0;
    std::vector<double> polarity_magnitudes;
    std::vector<double> polarity_angles;
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
            mean_magnitude += cell_iter->GetCellData()->GetItem("Polarity magnitude");
            mean_angle += cell_iter->GetCellData()->GetItem("Polarity angle");
            polarity_magnitudes.push_back(cell_iter->GetCellData()->GetItem("Polarity magnitude"));
            double angle = cell_iter->GetCellData()->GetItem("Polarity angle");
            if (angle>=EcadTurnoverUtilities::PI*3/4&&angle<=EcadTurnoverUtilities::PI*7/4)
            {
                angle +=EcadTurnoverUtilities::PI;
                if (angle>2*EcadTurnoverUtilities::PI)
                {
                    angle = angle - 2*EcadTurnoverUtilities::PI;
                }
            }
            polarity_angles.push_back(angle);


            num_cells++;
        }
    }

    double sum_sin = 0;
    double sum_cos = 0;
    for (unsigned int i=0; i<num_cells; ++i)
    {
        sum_sin += std::sin(polarity_angles[i]);
        sum_cos += std::cos(polarity_angles[i]);
    }
    mean_angle = std::atan2(sum_sin/num_cells, sum_cos/num_cells);

    const double average_x_sq = std::pow(sum_cos/num_cells,2);
    const double average_y_sq = std::pow(sum_sin/num_cells,2);
    const double mean_point_mag = std::sqrt(average_x_sq + average_y_sq);
    const double standard_dev = std::sqrt(-2*std::log(mean_point_mag));
    const double angle_sem = standard_dev/std::sqrt((double)num_cells);

    mean_magnitude /=num_cells;
    //mean_angle /= num_cells;
    double sum_squares_magnitude= 0, sum_squares_angle=0;
    for (unsigned int i=0; i<num_cells; ++i)
    {
        sum_squares_magnitude+=std::pow(polarity_magnitudes[i]-mean_magnitude,2);
        //sum_squares_angle+=std::pow(polarity_angles[i]-mean_angle,2);
    }
    double error_magnitudes = std::sqrt(sum_squares_magnitude/(num_cells-1))/std::sqrt((double)num_cells);
    //double error_angles = std::sqrt(sum_squares_angle/(num_cells-1))/std::sqrt((double)num_cells);
    *this->mpOutStream << mean_magnitude <<"\t";
    *this->mpOutStream << error_magnitudes <<"\t";
    *this->mpOutStream << mean_angle <<"\t";
    *this->mpOutStream << angle_sem <<"\t";
}

// Explicit instantiation
template class PolarityWriter<1,1>;
template class PolarityWriter<1,2>;
template class PolarityWriter<2,2>;
template class PolarityWriter<1,3>;
template class PolarityWriter<2,3>;
template class PolarityWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PolarityWriter)
