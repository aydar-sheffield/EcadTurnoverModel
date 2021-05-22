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

#include "EcadJunctionSrnModel.hpp"
#include "EcadJunctionOdeSystem.hpp"
#include "EcadTurnoverUtilities.hpp"
EcadJunctionSrnModel::EcadJunctionSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(3, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<EcadJunctionSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<EcadJunctionOdeSystem, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

EcadJunctionSrnModel::EcadJunctionSrnModel(const EcadJunctionSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
    assert(rModel.GetOdeSystem());
    AbstractOdeSystem* p_parent_system(rModel.GetOdeSystem());
    SetOdeSystem(new EcadJunctionOdeSystem(p_parent_system->rGetStateVariables()));
    for (unsigned int i=0; i < p_parent_system->GetNumberOfParameters(); ++i)
        mpOdeSystem->SetParameter(i, p_parent_system->GetParameter(i));


    /*mpOdeSystem->SetStateVariable("Mobile Ecad",EcadTurnoverUtilities::mobile_eq);
    mpOdeSystem->SetStateVariable("Bound Ecad",EcadTurnoverUtilities::bound_eq*0.3603);
    mpOdeSystem->SetStateVariable("Rest length",-1.0);*/
}

AbstractSrnModel* EcadJunctionSrnModel::CreateSrnModel()
{
    return new EcadJunctionSrnModel(*this);
}

void EcadJunctionSrnModel::SimulateToCurrentTime()
{
    // Custom behaviour
    UpdateKinetics();
    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void EcadJunctionSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new EcadJunctionOdeSystem);
}

void EcadJunctionSrnModel::InitialiseDaughterCell()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);
    mpOdeSystem->SetStateVariable("Mobile Ecad",EcadTurnoverUtilities::mobile_eq);
    mpOdeSystem->SetStateVariable("Bound Ecad",EcadTurnoverUtilities::bound_eq*0.36);
    mpOdeSystem->SetStateVariable("Rest length",-1.0);

    mpOdeSystem->SetParameter("Endocytosis rate",0.0);
    mpOdeSystem->SetParameter("Current length",0.015);
    mpOdeSystem->SetParameter("Remodelling rate",0.05);
}

void EcadJunctionSrnModel::UpdateKinetics()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double endo_rate
    = mpCell->GetCellEdgeData()->GetItem("Endocytosis rate")[this->GetEdgeLocalIndex()];
    double current_length
    = mpCell->GetCellEdgeData()->GetItem("Current length")[this->GetEdgeLocalIndex()];
    double k_L
    = mpCell->GetCellEdgeData()->GetItem("Remodelling rate")[this->GetEdgeLocalIndex()];

    mpOdeSystem->SetParameter("Endocytosis rate", endo_rate);
    mpOdeSystem->SetParameter("Current length", current_length);
    mpOdeSystem->SetParameter("Remodelling rate", k_L);
}

double EcadJunctionSrnModel::GetMobileEcad()
{
    assert(mpOdeSystem != nullptr);
    const double val = mpOdeSystem->rGetStateVariables()[0];
    return val;
}

void EcadJunctionSrnModel::SetMobileEcad(const double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[0] = value;
}

double EcadJunctionSrnModel::GetBoundEcad()
{
    assert(mpOdeSystem != nullptr);
    const double val = mpOdeSystem->rGetStateVariables()[1];
    return val;
}

void EcadJunctionSrnModel::SetBoundEcad(const double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[1] = value;
}

double EcadJunctionSrnModel::GetRestLength()
{
    assert(mpOdeSystem != nullptr);
    const double val = mpOdeSystem->rGetStateVariables()[2];
    return val;
}

void EcadJunctionSrnModel::SetRestLength(const double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[2] = value;
}


void EcadJunctionSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

void EcadJunctionSrnModel::AddSrnQuantities(AbstractSrnModel *p_other_srn,
                                            const double scale)
{
    auto other_srn
    = static_cast<EcadJunctionSrnModel*>(p_other_srn);
    const double other_mobile = other_srn->GetMobileEcad();
    const double other_bound = other_srn->GetBoundEcad();
    const double this_mobile = GetMobileEcad();
    const double this_bound = GetBoundEcad();
    SetMobileEcad(0.5*(this_mobile+scale*other_mobile));
    SetBoundEcad(0.5*(this_bound+scale*other_bound));
}

void EcadJunctionSrnModel::AddShrunkEdgeSrn(AbstractSrnModel *p_shrunk_edge_srn)
{
    // Here we assume that one half of srn quantities are endocytosed and the remaining
    // half are split between neighbouring junctions. Hence we add 1/4 of srn variables
    // Note that this is not applicable to the rest length. Therefore, we store the
    // old value and copy it back
    //const double rest_length = GetRestLength();
    //AddSrnQuantities(p_shrunk_edge_srn, 0.00);
    //SetRestLength(rest_length);
}

void EcadJunctionSrnModel::AddMergedEdgeSrn(AbstractSrnModel* p_merged_edge_srn)
{
    // Add all srn variables to this edge srn
    // Note that this is not applicable to the rest length. Therefore, we store the
    // old value and copy it back
    //const double rest_length = GetRestLength();
    AddSrnQuantities(p_merged_edge_srn);
    SetRestLength(-1.0);
}

void EcadJunctionSrnModel::SplitEdgeSrn(const double relative_position)
{
    ScaleSrnVariables(1.0);
    const double current_rest_length = GetRestLength();
    SetRestLength(-1.0);
}


// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(EcadJunctionSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(EcadJunctionSrnModel)
