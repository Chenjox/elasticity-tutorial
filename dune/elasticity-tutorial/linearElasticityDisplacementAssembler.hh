#pragma once



#include <dune/istl/matrix.hh>
#include <dune/common/rangeutilities.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/assembler/denselocalassemblerbase.hh>

#include <dune/elasticity-tutorial/materials/material.hh>

using namespace Dune;

template <class LocalView, class Material, class RHSFunc>
class LocalLinearElasticityAssembler
  : public Assembler::DenseLocalAssemblerBase<LocalView>
{
    using GridView = typename LocalView::GridView;
    using Element = typename LocalView::Element;
    using LocalRHSFunc = typename RHSFunc::LocalFunction;

    static const int dim = Element::dimension;


private:
  RHSFunc rhs_;
  LocalRHSFunc localRhs_;
  Material material_;
  const LocalView* testLocalView_ = nullptr;
  const LocalView* trialLocalView_ = nullptr;
  std::vector<Dune::FieldMatrix<double,1,dim>> referenceJacobians_;
  std::vector<Dune::FieldMatrix<double,1,dim>> jacobians_;
  std::vector<Dune::FieldVector<double,1>> values_;
  const int quadRuleOrder_ = 2;

public: 

    template<class Basis>
    LocalLinearElasticityAssembler (const Basis&, const RHSFunc& rhs, const Material mat, const int quadRuleOrder)
        : rhs_(rhs),
        material_(mat),
        localRhs_(localFunction(rhs)),
        quadRuleOrder_(quadRuleOrder)
    {}

    void bindLocalViews (const LocalView& testLocalView, const LocalView& trialLocalView)
    {
      testLocalView_ = &testLocalView;
      trialLocalView_ = &trialLocalView;
    }

    void bindElement (const Element& element)
    {}

    template <class LocalMatrix>
    void assembleElementMatrix (const Element& element, LocalMatrix& localMatrix)
    {
      assert(testLocalView_ == trialLocalView_);

      // Get the Element Geometry
      auto geometry = element.geometry();
  
      // Get a Node, at which the Displacements happen.
      const auto& node = testLocalView_->tree();

      // Get a basis function.
      // Because we assume the node to be a power node of identical lagrange nodes
      // the first one will suffice.
      const auto& finiteElement = node.child(0).finiteElement();

      // Get the correct sizes for the jacobians
      // Probably this way to reduce allocs
      referenceJacobians_.resize(finiteElement.size());
      jacobians_.resize(finiteElement.size());
  
      // Get a Quadrature Rule
      const auto& quad = Dune::QuadratureRules<double, dim>::rule(element.type(), quadRuleOrder_);
  
      // Iteration over all Points in Questions
      for (const auto& [position, weight] : quad)
      {
        // Integration Factor because of variable substitution to the reference element.
        const auto integrationWeight = geometry.integrationElement(position) * weight;
  
        const auto geometryJacobianInverse = geometry.jacobianInverse(position);
        finiteElement.localBasis().evaluateJacobian(position, referenceJacobians_);

        // Derivatives in World Coordinates
        for (auto i : Dune::range(finiteElement.size()))
          jacobians_[i] = referenceJacobians_[i] * geometryJacobianInverse;

        // We will now do a Nodal Assembly Process
        //std::cout << jacobians_[0][0] << std::endl;
        // Real Space: ik, Test Space jl
        for (auto i : Dune::range(finiteElement.size()))
          for (auto k : Dune::range(dim))
            for (auto j : Dune::range(finiteElement.size()))
              for (auto l : Dune::range(dim)){
                // Get the correct Indizes of the local stiffness Matrix
                auto localRealIndex = node.child(k).localIndex(i);
                auto localVirtIndex = node.child(l).localIndex(j);

                auto displacementGradient = FieldMatrix<double, dim, dim>{0.0};
                // Make a Displacement Gradient
                for (auto ii : Dune::range(dim))
                  displacementGradient[k][ii] = jacobians_[i][0][ii];
                // Get the actual Strains
                auto strains = FieldMatrix<double, dim, dim>{0.0};

                for (auto ii : Dune::range(dim))
                  for (auto jj : Dune::range(dim))
                    strains[ii][jj] = 0.5*(displacementGradient[ii][jj] + displacementGradient[jj][ii]);

                auto stresses = material_.stresses(strains);
                
                // Same for the virtual strain
                displacementGradient = 0.0;

                for (auto ii : Dune::range(dim))
                  displacementGradient[l][ii] = jacobians_[j][0][ii];
                auto virtStrains = FieldMatrix<double, dim, dim>{0.0};

                for (auto ii : Dune::range(dim))
                  for (auto jj : Dune::range(dim))
                    virtStrains[ii][jj] = 0.5*(displacementGradient[ii][jj] + displacementGradient[jj][ii]);

                // Contraction
                auto value = 0.0;
                for (auto ii : Dune::range(dim))
                  for (auto jj : Dune::range(dim))
                    value += stresses[ii][jj] * virtStrains[ii][jj];

                localMatrix[localRealIndex][localVirtIndex] += value * integrationWeight;
            
              }
      }
    }

     void bindLocalView (const LocalView& testLocalView)
    {
      testLocalView_ = &testLocalView;
    }


    template <class LocalVector>
  void assembleElementVector (const Element& element, LocalVector& localVector)
  {
    auto geometry = element.geometry();

    localRhs_.bind(element);

    const auto& node = testLocalView_->tree();
    const auto& finiteElement = node.child(0).finiteElement();
    values_.resize(finiteElement.size());

    const auto& quad = Dune::QuadratureRules<double, dim>::rule(element.type(), quadRuleOrder_);

    for (const auto& [position, weight] : quad)
    {
      const auto integrationWeight = geometry.integrationElement(position) * weight;

      finiteElement.localBasis().evaluateFunction(position, values_);

      auto rhsValue = localRhs_(position);

      for (auto k : Dune::range(dim))
        for (auto i : Dune::range(finiteElement.size()))
          localVector[node.child(k).localIndex(i)] += values_[i] * rhsValue[k] * integrationWeight;
    }
  }


};


template <class Basis, class RHSFunc, class Material>
LocalLinearElasticityAssembler(const Basis&, const RHSFunc& rhs, const Material mat, const int quadRuleOrder) -> LocalLinearElasticityAssembler<typename Basis::LocalView, Material, RHSFunc>;
