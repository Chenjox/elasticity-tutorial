// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "dune/common/bitsetvector.hh"
#include "dune/common/fvector.hh"
#include "dune/elasticity-tutorial/materials/stvenantkirchhoff.hh"
#include "dune/functions/gridfunctions/gridviewfunction.hh"
#include "dune/istl/bvector.hh"
#include <cmath>
#include <ostream>
#include <string>
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/parametertree.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/version.hh>
#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/uggrid.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>

#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>

#include <dune/assembler/defaultglobalassembler.hh>
#include <dune/assembler/backends/istlvectorbackend.hh>
#include <dune/assembler/backends/istlmatrixbackend.hh>

#include <dune/fufem/dunepython.hh>
#include <dune/gmsh4/gmsh4reader.hh>

#include <dune/solvers/iterationsteps/cgstep.hh>
#include <dune/solvers/iterationsteps/istlseqilu0step.hh>


#include <dune/solvers/solvers/loopsolver.hh>
#include <dune/solvers/norms/twonorm.hh>

#include <dune/elasticity-tutorial/linearElasticityDisplacementAssembler.hh>
#include <dune/elasticity-tutorial/statistics.hh>
#include <dune/elasticity-tutorial/meshUtilities.hh>

#include <dune/vtk/vtkwriter.hh>
/////////////
// Loop over LagrangeNodes
/////////////

// grid dimension
const int dim = 2;
const int dimworld = 2;

const int displacementOrder = 2;

using namespace Dune;

int main(int argc, char** argv)
{
  // Maybe initialize MPI
  Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);

  ////////////////////////////////////////////////
  // Parameter File Reading and Sanity Checking
  ////////////////////////////////////////////////

  // Check for appropriate number of command line arguments
  if (argc < 3)
    DUNE_THROW(Exception, "Usage: ./elasticity-tutorial <local path to cases> <case python file without .py>");


  // Start Python interpreter

  Python::start();
  auto pyMain = Python::main();
  Python::run("import math");

  //feenableexcept(FE_INVALID);

  Python::runStream()
    << std::endl << "import sys"
    << std::endl << "import inspect"
    << std::endl << "from pathlib import Path"
    << std::endl << "sys.path.append('" << argv[1] << "')"
    << std::endl;

  // parse data file
  auto pyModule = pyMain.import(argv[2]);

  // Get main parameter set
  ParameterTree parameterSet;
  pyModule.get("parameterSet").toC(parameterSet);

  ParameterTreeParser::readOptions(argc, argv, parameterSet);

  const double firstLameParameter = parameterSet.get<double>("firstLameParameter");
  const double secondLameParameter= parameterSet.get<double>("secondLameParameter");
  const int refinementLevels = parameterSet.get<int>("refinementLevels",0);
  const int quadRuleOrder = parameterSet.get<int>("quadRuleOrder",2);
  const std::string filename = parameterSet.get<std::string>("gridFile");

  //////////////////////////////////
  // Grid Setup
  //////////////////////////////////

  std::cout << "Reading Grid file: " << filename << std::endl;

  using Grid = UGGrid<dimworld>;
  std::unique_ptr grid = Gmsh4Reader<Grid>::createGridFromFile(filename);
  grid->globalRefine(refinementLevels);

  using GridView = typename Grid::LeafGridView;
  auto gridView = grid->leafGridView();

  Tutorial::checkMesh(gridView);

  //////////////////////////////////
  // Choosing an Ansatz Space
  //////////////////////////////////

  using namespace Functions::BasisFactory;

  auto displacementBasis = makeBasis(
    gridView,
    power<dim>(
      lagrange<displacementOrder>()
    ));

  
  using LocalView = typename decltype(displacementBasis)::LocalView;

  /////////////////////////////////////////////////////////
  // Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////

  using Vector = BlockVector<FieldVector<double,dim>>;
  using Matrix = BCRSMatrix<FieldMatrix<double,dim,dim>>;

  // Indicator sets (= std::array<char,dim> für inner und std::vector für outer)
  using BitVector = BitSetVector<dim>;

  /////////////////////////////////////////////////////////
  // Assemble the system
  /////////////////////////////////////////////////////////

  Vector rhs;
  auto rhsBackend = Assembler::ISTLVectorBackend(rhs);
  //Functions::istlVectorBackend(rhs);
  rhsBackend.resize(displacementBasis);
  rhs = 0;
  Matrix stiffnessMatrix;

  /////////////////////////////////////////////////////////
  // Hack to get the true position of Lagrange Nodes
  /////////////////////////////////////////////////////////

  Vector positionOfNodes;
  auto positionOfNodesBackend = Functions::istlVectorBackend(positionOfNodes);
  positionOfNodesBackend.resize(displacementBasis);

  Functions::interpolate(displacementBasis, positionOfNodesBackend, [&] (auto&& xi) { return xi; });

  //////////////////////////////////////////////////////////

  // Hier wird die Funktion in der Datei zu einer C++ Funktion gemacht. Finde ich besser als ein lambda.
  // EDIT: Soll wohl sehr langsam sein.
  auto dirichletIndicator = Python::make_function<FieldVector<bool,dim>>(pyModule.get("DirichletIndicatorFunction"));

  BitVector isDirichletBoundary;
  auto isDirichletBoundaryBackend = Functions::istlVectorBackend(isDirichletBoundary);
  isDirichletBoundaryBackend.resize(displacementBasis);

  // Initialisierung
  for (auto&& b0i : isDirichletBoundary)
    for (std::size_t j=0; j<b0i.size(); ++j)
      b0i[j] = false;

  // Mark each Boundary Intersection in the Vector.
  // The Callback gets a globalIndex corresponding to a degree of freedom
  // Therefore it is sufficient to simply Index the positionOfNodes Vector for its positions
  // as the Indices will be associated to the nodes.
  Functions::forEachBoundaryDOF(displacementBasis, [&] (auto && globalIndex){
    // Each Node is indicated by a function depending on the position of the node.
    // It returns a vector of bool, therefore it needs to be indexed by the second index.
    isDirichletBoundaryBackend[globalIndex] = dirichletIndicator(positionOfNodes[globalIndex[0]])[globalIndex[1]];
  });
/* Das geht eventuelle auch....
auto predicate = [](auto x)
{
return x[0] < 1e-8
|| x[1] < 1e-8
|| (x[0] > 0.4999 && x[1] > 0.4999);
};
// Evaluating the predicate will mark all Dirichlet degrees of freedom
std::vector<bool> dirichletNodes;
Functions::interpolate(basis, dirichletNodes, predicate);
*/
  // Dirichlet boundary value function, depending on the homotopy parameter
  // It is a Python Class, that gets instantiation with a parameter
  Python::Callable dirichletBoundaryPythonClass = pyModule.get("DirichletBoundary");
  // Construct with a particular value of the homotopy parameter
  Python::Reference dirichletValuesPythonObject = dirichletBoundaryPythonClass(0.0);

  
  auto dirichletValueFunction = Python::make_function<FieldVector<double,dim>>(dirichletValuesPythonObject.get("dirichletValues"));

  // Create the actual Dirichlet values in the lefthandside
  Functions::interpolate(
    displacementBasis,
    rhs,
    dirichletValueFunction,
    isDirichletBoundary
  );
  
  // The force in Question.
  auto bodyforceFunction = Python::make_function<FieldVector<double,dim>>(pyModule.get("bodyforce"));

  ////////////////////////////////////////////////
  // Assembly
  ////

  auto rhsFunction = Functions::makeGridViewFunction(bodyforceFunction,gridView);

  // Ein Backend, daraus eine Pattern
  auto stiffnessMatrixBackend = Assembler::ISTLMatrixBackend(stiffnessMatrix);
  auto stiffnessMatrixPattern = stiffnessMatrixBackend.patternBuilder();
  stiffnessMatrixPattern.resize(displacementBasis,displacementBasis);
  

  // Global Assembler
  auto assembler = Assembler::Assembler(displacementBasis);

  // Local Assember
  auto mat = Tutorial::StVenantKirchhoff<dimworld>{firstLameParameter, secondLameParameter};
  //
  auto localAssembler = LocalLinearElasticityAssembler{ displacementBasis, rhsFunction, mat, quadRuleOrder };

  assembler.assembleMatrixPattern(localAssembler, stiffnessMatrixPattern);
  stiffnessMatrixPattern.setupMatrix();

  stiffnessMatrix = 0;

  assembler.assembleMatrix(localAssembler, stiffnessMatrixBackend);
  assembler.assembleVector(localAssembler, rhsBackend);

  //////////////////////////////////////////////////
  // Solution Algorithm.
  /////////////////////////////////////////////////

  double stepTol    = 1e-12; // termination criterion
  double solveTol   = 1e-10;  // error in the solution
  auto norm = TwoNorm<Vector>();

  Vector x;
  auto xBackend = Functions::istlVectorBackend(x);
  xBackend.resize(displacementBasis);
  xBackend = 0;

  //
  //auto precon = ISTLSeqILU0Step<Matrix, Vector>(1.0);

  // No Preconditioner
  Solvers::CGStep<Matrix,Vector> cgStep(stiffnessMatrix, x, rhs);
  cgStep.setIgnore(isDirichletBoundary);

  Solvers::LoopSolver<Vector> solver(cgStep, 500, stepTol, norm, Solver::QUIET, false);

  solver.check();
  solver.preprocess();
  solver.solve();

  auto solutionFunction = Functions::makeDiscreteGlobalBasisFunction<FieldVector<double,dim>>(displacementBasis,x);

  Vtk::VtkWriter<GridView> vtkWriter(gridView, Vtk::FormatTypes::ASCII, Vtk::DataTypes::FLOAT64);
  vtkWriter.addPointData(solutionFunction, "displacement");

  std::string path = argv[1];
  std::string t = argv[2];
  vtkWriter.write(path + "/" + t+".vtu");
}
