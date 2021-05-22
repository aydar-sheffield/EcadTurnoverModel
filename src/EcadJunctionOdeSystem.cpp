/*

Copyright (c) 2005-2020, University of Oxford.
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

#include "CellwiseOdeSystemInformation.hpp"
#include "EcadJunctionOdeSystem.hpp"
#include "EcadTurnoverUtilities.hpp"

EcadJunctionOdeSystem::EcadJunctionOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(3)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<EcadJunctionOdeSystem>);

    SetDefaultInitialCondition(0, 1.0);
    SetDefaultInitialCondition(1, 1.0);
    SetDefaultInitialCondition(2, 0.05);
    //By default zero. If no interior SRN model is specified, interior parameters is zero
    for (unsigned int i=0; i<3; ++i)
        this->mParameters.push_back(0.0);
    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

EcadJunctionOdeSystem::~EcadJunctionOdeSystem()
{
}

void EcadJunctionOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    const double mobile = rY[0];
    const double bound = rY[1];
    const double k_minus = this->mParameters[0];

    const double binding = EcadTurnoverUtilities::k_plus*mobile;
    const double recruitment = EcadTurnoverUtilities::k_c_plus;
    const double mobile_degradation = EcadTurnoverUtilities::k_c_minus*mobile;
    const double endocytosis = k_minus*bound;

    rDY[0] = recruitment - mobile_degradation - binding;
    rDY[1] = binding - endocytosis;

    //Rest length
    const double rest_length = rY[2];
    const double current_length = this->mParameters[1];
    const double k_L = this->mParameters[2];
    rDY[2] = k_L*(current_length-rest_length);
}

template<>
void CellwiseOdeSystemInformation<EcadJunctionOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Mobile Ecad");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Bound Ecad");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Rest length");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.05);

    this->mParameterNames.push_back("Endocytosis rate");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Current length");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("Remodelling rate");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(EcadJunctionOdeSystem)
