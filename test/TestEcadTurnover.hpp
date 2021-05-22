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

#ifndef TESTT1AFTERGROWTH_HPP_
#define TESTT1AFTERGROWTH_HPP_


#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"

#include <chrono>

#include "EcadTurnoverUtilities.hpp"
#include "EcadTurnoverForce.hpp"
#include "RandomG1CellCycle.hpp"
#include "PreStretchModifier.hpp"
#include "PostStretchModifier.hpp"
#include "EcadTurnoverModifier.hpp"
#include "EcadJunctionSrnModel.hpp"
#include "AverageElongationWriter.hpp"
#include "AverageAreaWriter.hpp"
#include "AverageTensionWriter.hpp"
#include "T1SwapWriter.hpp"
#include "EcadWriter.hpp"
#include "VertexMeshReader.hpp"


#include "VtkVertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"

#include "PolygonalDistributionWriter.hpp"
#include "JunctionalLengthWriter.hpp"
#include "AreaAndBoundaryWriter.hpp"
#include "PolarityWriter.hpp"
#include "AverageEcadWriter.hpp"

#include "CommandLineArguments.hpp"
/* Having included all the necessary header files, we proceed by defining the test class.
 */
class TestEcadTurnover : public AbstractCellBasedTestSuite
{
public:

    void TestMechanoSensitiveEcad()
    {
        RandomNumberGenerator::Instance()->Reseed(time(NULL));
        /* We include the next line because Vertex simulations cannot be run in parallel */
        EXIT_IF_PARALLEL;
        unsigned int outp = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-opt");
        std::cout<<"Stretching with mesh # "<<outp<<std::endl;
        //VertexMeshReader<2,2> mesh_reader("../../Test_output/TestRandomMeshGeneration/vertex_mesh_800"); prior to automatisation
        VertexMeshReader<2,2> mesh_reader("../../EcadSimulations/Meshes/vertex_mesh_800_"+std::to_string(outp));
        MutableVertexMesh<2,2> vertex_mesh;
        vertex_mesh.ConstructFromMeshReader(mesh_reader);
        MutableVertexMesh<2,2>* p_mesh(&vertex_mesh);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_diff_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_nondividing_type);
        for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            /* Initialise cell cycle */
            RandomG1CellCycle* p_cc_model = new RandomG1CellCycle();
            p_cc_model->SetMaxTransitGenerations(0);
            p_cc_model->SetTransitCellG1Duration(EcadTurnoverUtilities::G1AverageDuration);
            p_cc_model->SetG2Duration(EcadTurnoverUtilities::G2Duration);
            p_cc_model->SetSDuration(1e-12);
            p_cc_model->SetMDuration(1e-12);
            p_cc_model->SetDimension(2);

            auto p_cell_srn_model = new CellSrnModel();
            for (unsigned i = 0; i < p_mesh->GetElement(elem_index)->GetNumEdges(); i ++)
            {
                MAKE_PTR(EcadJunctionSrnModel, p_edge_model);
                double mobile = EcadTurnoverUtilities::mobile_eq;
                double bound = EcadTurnoverUtilities::bound_eq*0.3603;
                double current_length = p_mesh->GetElement(elem_index)->GetEdge(i)->rGetLength();
                std::vector<double> init_conditions = {mobile, bound, current_length};

                p_edge_model->SetInitialConditions(init_conditions);
                p_cell_srn_model->AddEdgeSrnModel(p_edge_model);
            }

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_srn_model));
            p_cell->SetCellProliferativeType(p_nondividing_type);

            double birth_time
            = -RandomNumberGenerator::Instance()->ranf()*(EcadTurnoverUtilities::G1AverageDuration+EcadTurnoverUtilities::G2Duration);
            birth_time *= 0.25;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.InitialiseCells();
        cell_population.SetDampingConstantNormal(EcadTurnoverUtilities::friction);

        cell_population.AddPopulationWriter<AverageElongationWriter>();
        cell_population.AddPopulationWriter<AverageAreaWriter>();
        cell_population.AddPopulationWriter<AverageTensionWriter>();
        cell_population.AddPopulationWriter<EcadWriter>();
        cell_population.AddPopulationWriter<AreaAndBoundaryWriter>();
        cell_population.AddPopulationWriter<PolarityWriter>();
        cell_population.AddPopulationWriter<AverageEcadWriter>();

        // Some cell data have to be initialised prior to stretching
        // For example, we need to know which cells are to be stretched
        MAKE_PTR(PreStretchModifier<2>, p_pre_modifier);
        p_pre_modifier->UpdateCellData(cell_population);

        MAKE_PTR(EcadTurnoverModifier<2>, p_ecad_modifier);
        p_ecad_modifier->UpdateCellData(cell_population);

        OffLatticeSimulation<2> simulator(cell_population, false, false);
        simulator.SetOutputDirectory("TestEcadTurnover/sim_number_"+std::to_string(outp));

        //This is a topology update modifier, because rest length values after
        // remeshing will need to be updated so that vertex time-stepping is updated properly
        simulator.AddTopologyUpdateSimulationModifier(p_ecad_modifier);

        MAKE_PTR(EcadTurnoverForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(TargetAreaLinearGrowthModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        cell_population.SetWriteCellVtkResults(false);
        cell_population.SetWriteEdgeVtkResults(false);


        try
        {
            std::clock_t c_start = std::clock();
            auto t_start = std::chrono::high_resolution_clock::now();

            //Initial 5 minute stretch requires smaller step size
            simulator.SetSamplingTimestepMultiple(800);
            simulator.SetEndTime(EcadTurnoverUtilities::stretch_time);
            simulator.SetDt(EcadTurnoverUtilities::dt/8);
            simulator.Solve();

            //After the tissue had been stretched
            simulator.SetSamplingTimestepMultiple(100);
            simulator.SetEndTime(EcadTurnoverUtilities::Tend);
            simulator.SetDt(EcadTurnoverUtilities::dt);
            simulator.Solve();

            std::clock_t c_end = std::clock();
            auto t_end = std::chrono::high_resolution_clock::now();
            const double cpu_time = (c_end-c_start)/CLOCKS_PER_SEC;
            const double wall_time = std::chrono::duration_cast<std::chrono::minutes>(t_end - t_start).count();
            std::cout<<"CPU time: "<<cpu_time/60.0<<" minutes"<<std::endl;
            std::cout<<"Wall time: "<<wall_time<<" minutes"<<std::endl;
        }
        catch(const std::exception& exc)
        {
            std::cerr<<exc.what()<<std::endl;
        }
    }


};


#endif
