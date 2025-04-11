#pragma once

#include "dune/elasticity-tutorial/materials/material.hh"
#include <dune/common/fmatrix.hh>

namespace Dune::Tutorial {


template<const int dim>
class StVenantKirchhoff : Material<dim> {
    private:
        double firstLameParameter_, secondLameParameter_;

    public:
    
    StVenantKirchhoff(double firstLameParameter, double secondLameParameter):
        firstLameParameter_(firstLameParameter),
        secondLameParameter_(secondLameParameter) {};

    FieldMatrix<double,dim,dim> stresses(FieldMatrix<double,dim,dim> strains){

        FieldMatrix<double, dim, dim> stresses{0.0};

        double trace = 0.0;
        for (int i : Dune::range(dim)){
            trace += strains[i][i];
            for (int j : Dune::range(dim)){
                stresses[i][j] = 2.0 * secondLameParameter_ * strains[i][j];
            }
        }
        for (int i : Dune::range(dim)){
            stresses[i][i] += firstLameParameter_ * (trace / 3.0);
        }

        return stresses;
    };

};

}