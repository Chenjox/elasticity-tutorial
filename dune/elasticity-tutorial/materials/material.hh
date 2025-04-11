#pragma once

#include <dune/common/fmatrix.hh>

namespace Dune::Tutorial {


template<const int dim>
class Material {

    virtual FieldMatrix<double,dim,dim> stresses(FieldMatrix<double,dim,dim> strains);

};

}