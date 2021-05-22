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
#include "FarhadifarForce.hpp"

#include <chrono>

#include "EcadTurnoverUtilities.hpp"
#include "RandomG1CellCycle.hpp"
#include "VertexMeshWriter.hpp"

#include "CommandLineArguments.hpp"
/**
 * This test generates a mesh of approximately 800 cells
 */
class TestGenerateRandomEpithelia : public AbstractCellBasedTestSuite
{
public:

    void TestGenerateRandomMesh()
    {
        RandomNumberGenerator::Instance()->Reseed(time(NULL));
        /* We include the next line because Vertex simulations cannot be run in parallel */
        EXIT_IF_PARALLEL;
        HoneycombVertexMeshGenerator generator(5, 5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        unsigned int outp = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-opt");
        std::cout<<"Generating mesh # "<<outp<<std::endl;

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_diff_type);
        for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
             /*Initialise cell cycle*/
            RandomG1CellCycle* p_cc_model = new RandomG1CellCycle();
            p_cc_model->SetMaxTransitGenerations(4);
            p_cc_model->SetTransitCellG1Duration(EcadTurnoverUtilities::G1AverageDuration);
            p_cc_model->SetG2Duration(EcadTurnoverUtilities::G2Duration);
            p_cc_model->SetSDuration(1e-12);
            p_cc_model->SetMDuration(1e-12);
            p_cc_model->SetDimension(2);


            CellPtr p_cell(new Cell(p_state, p_cc_model));
            p_cell->SetCellProliferativeType(p_diff_type);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.SetDampingConstantNormal(EcadTurnoverUtilities::friction);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestEcadTurnoverMesh");
        simulator.SetSamplingTimestepMultiple(500);

        simulator.SetEndTime(50*EcadTurnoverUtilities::real_to_simulated);
        simulator.SetDt(0.01);

        MAKE_PTR(FarhadifarForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(TargetAreaLinearGrowthModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        const std::string mesh_name = "vertex_mesh_800_"+std::to_string(outp);
        try
        {
            std::clock_t c_start = std::clock();
            auto t_start = std::chrono::high_resolution_clock::now();

            simulator.Solve();
            // Create a vertex mesh writer
            VertexMeshWriter<2,2> vertex_mesh_writer("TestRandomMeshGeneration", mesh_name, false);

            vertex_mesh_writer.WriteFilesUsingMesh(*p_mesh);

            std::clock_t c_end = std::clock();
            auto t_end = std::chrono::high_resolution_clock::now();
            const double cpu_time = (c_end-c_start)/CLOCKS_PER_SEC;
            const double wall_time = std::chrono::duration_cast<std::chrono::minutes>(t_end - t_start).count();
            std::cout<<"CPU time: "<<cpu_time/60.0<<" minutes"<<std::endl;
            std::cout<<"Wall time: "<<wall_time<<" minutes"<<std::endl;
            std::cout<<"Total number of cells: "<<p_mesh->GetNumElements()<<std::endl;
        }
        catch(const std::exception& exc)
        {
            std::cerr<<exc.what()<<std::endl;
        }
    }


};


#endif
