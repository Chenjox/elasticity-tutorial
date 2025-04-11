#pragma once


#include <algorithm>
#include <cmath>
#include <limits>


namespace Dune::Tutorial
{


  // https://math.stackexchange.com/questions/2148877/iterative-calculation-of-mean-and-standard-deviation
    class StatisticsAccumulator
    {
        double _minimum; 
        double _maximum;
        double _sumOfValues;
        double _sumOfSquaredValues;
        double _count;

        public:

            StatisticsAccumulator(){
                _minimum = std::numeric_limits<double>::max(); // Alles sollte kleiner als das Maximum sein
                _maximum = std::numeric_limits<double>::min(); // Alles größer als dieses Minimum
                _sumOfValues = 0.0;
                _sumOfSquaredValues = 0.0;
                _count = 0.0;
            }

            void add_datum(double value){
                _minimum = std::min(_minimum, value);
                _maximum = std::max(_maximum, value);
                _sumOfValues += value;
                _sumOfSquaredValues += value*value;
                _count += 1.0;
            }

            double calculate_mean(){
                return _sumOfValues/_count;
            }

            double calculate_std_deviation(){
                return std::sqrt( _sumOfSquaredValues/_count - (_sumOfValues/_count)*(_sumOfValues/_count) );
            }

            double calculate_min(){
                return _minimum;
            }

            double calculate_max(){
                return _maximum;
            }
    };

}