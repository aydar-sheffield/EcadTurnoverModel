/*
 * RandomG1CellCycle.cpp
 *
 *  Created on: 21 Nov 2019
 *      Author: aydar
 */

#include "RandomG1CellCycle.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

RandomG1CellCycle::RandomG1CellCycle()
{
}

RandomG1CellCycle::RandomG1CellCycle(const RandomG1CellCycle& rModel)
   : AbstractSimpleGenerationalCellCycleModel(rModel)
{
    /*
     * The member variables mGeneration and mMaxTransitGeneration are
     * initialized in the AbstractSimpleGenerationalCellCycleModel
     * constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* RandomG1CellCycle::CreateCellCycleModel()
{
    return new RandomG1CellCycle(*this);
}

void RandomG1CellCycle::SetG1Duration()
{
    /**
     * Here we focus only on transit cell proliferative type. Other types use default parameter values
     */
    assert(mpCell != nullptr);
    double u = RandomNumberGenerator::Instance()->ranf();
    if ( mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>() )
    {
        mG1Duration = GetStemCellG1Duration() + 4*u; // U[14,18] for default parameters (mStemCellG1Duration) according to Meineke
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        mG1Duration = -log(u)*GetTransitCellG1Duration();
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

void RandomG1CellCycle::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractSimpleGenerationalCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

