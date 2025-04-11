#pragma once

#include <algorithm>


#include <iostream>
#include <limits>
#include <vector>

#include <dune/elasticity-tutorial/statistics.hh>

namespace Dune::Tutorial 
{

    template<class Geometry>
    static auto diameter (const Geometry& geometry)
    {
      typename Geometry::ctype h = 0.0;
      for (int i=0; i<geometry.corners(); i++)
        for (int j=i+1; j<geometry.corners(); j++)
          h = std::max(h, (geometry.corner(j)-geometry.corner(i)).two_norm());
      return h;
    }



    template<class Geometry>
    static auto nearestPointProjectionOntoConvexHull(const Geometry& geom, std::vector<int> &convexHull){
        typename Geometry::GlobalCoordinate iterate(0);
        for (int i=0; i< (int) convexHull.size(); i++)
            iterate += geom.corner(convexHull[i]);
        iterate /= convexHull.size();

        int k = 0;

        while (true && k < 1000) {
            std::vector<int> maxDistPointIndezes{};
            maxDistPointIndezes.reserve(convexHull.size()); // die Convexe Hülle von z kann maximal z elemente besitzen.
            typename Geometry::ctype maxDist = 0.0;

            // Bestimmen der Menge an Punkten mit maximalem Abstand
            for (int i=0; i< (int) convexHull.size(); i++){
                double dist = (iterate - geom.corner(convexHull[i])).two_norm();
                if (dist > maxDist || std::abs(dist - maxDist) <= 1e-14){
                    maxDistPointIndezes.push_back(i);
                    maxDist = dist;
                }
            }
            // Sind elemente der Convexen Hülle gleich mit der Abstandsmenge?
            if (convexHull.size() == maxDistPointIndezes.size()){
                bool check = true; // Nehme an, es sei gleich
                // Da beide Listen in gleicher weise Befüllt werden, ist die Ordnung der Indizes der Listen gleich,
                // bei gleicher größe müssen also die Indizes auch in gleicher Reihenfolge vorkommen
                for (int i=0; i < (int) convexHull.size(); i++){
                    check = check && convexHull[i] == maxDistPointIndezes[i];
                }
                if (check) {
                    return iterate;
                }
            }
            // Wenn dem nicht so ist, dann REKURSION!
            auto yk = nearestPointProjectionOntoConvexHull(geom, maxDistPointIndezes);
            // Anschließend wird alpha k bestimmt.
            double alpha = std::numeric_limits<double>::max();
            for (int i = 0; i < (int) convexHull.size(); i++){
                // Wenn der Eintrag an der Stelle i **nicht** enthalten ist
                if (std::find(maxDistPointIndezes.begin(), maxDistPointIndezes.end(), convexHull[i]) == maxDistPointIndezes.end()){
                    double under =(yk - iterate).dot( geom.corner(convexHull[i]) - yk );
                    if (under < 0.0){
                        continue;
                    }
                    alpha = std::min(alpha, (geom.corner(convexHull[i]) - iterate).two_norm2() - maxDist*maxDist)
                    / ( 2.0 * under );
                }
            }
            if (alpha >= 1.0){
                return iterate;
            }else {
                iterate = iterate + alpha*(yk - iterate);
            }

            k+=1;
        }
        // this should be unreachable
        return iterate;
    }


    template<class Geometry>
    static auto inscribedRadius(const Geometry& geometry)
    {
        //typename Geometry::ctype inRadius = 0.0;
        // Outer sphere
        // Algorithm from "An Algorithm for Finding the Chebyshev Center of a Convex Polyhedron" von N. D. Botkin and V. L. Turova-Botkina
        typename Geometry::GlobalCoordinate iterate(0);
        for (int i=0; i<geometry.corners(); i++)
            iterate += geometry.corner(i);
        iterate /= geometry.corners();

        // Buchführung
        typename Geometry::ctype maxDist = 0.0;
        int k = 0;
        while (true && k < 1000 ) {
            // Schritt 1: Finde Punkt mit Maximalem Abstand in den Punkten des Elements
            std::vector<int> maxDistPointIndezes{};
            maxDistPointIndezes.reserve(geometry.corners()); // die Convexe Hülle von z kann maximal z elemente besitzen.

            // Bestimmen der Menge an Punkten mit maximalem Abstand
            for (int i=0; i<geometry.corners(); i++){
                double dist = (iterate - geometry.corner(i)).two_norm();
                if (dist > maxDist || std::abs(dist - maxDist) <= 1e-14){
                    maxDistPointIndezes.push_back(i);
                    maxDist = dist;
                }
            }

            if (geometry.corners() == (int) maxDistPointIndezes.size()){
                bool check = true; // Nehme an, es sei gleich
                // Da beide Listen in gleicher weise Befüllt werden, ist die Ordnung der Indizes der Listen gleich,
                // bei gleicher größe müssen also die Indizes auch in gleicher Reihenfolge vorkommen
                for (int i=0; i < geometry.corners(); i++){
                    check = check && i == maxDistPointIndezes[i];
                }
                if (check) {
                    break;
                    // return iterate;
                }
            }
            // Wenn dem nicht so ist, dann REKURSION!
            auto yk = nearestPointProjectionOntoConvexHull(geometry, maxDistPointIndezes);

            // Nächste Punkt Projektion auf die Convexe Hülle von maxDistPointIndezes.
            double alpha = std::numeric_limits<double>::max();
            for (int i = 0; i < (int) geometry.corners(); i++){
                // Wenn der Eintrag an der Stelle i **nicht** enthalten ist
                if (std::find(maxDistPointIndezes.begin(), maxDistPointIndezes.end(), i) == maxDistPointIndezes.end()){
                    double under =(yk - iterate).dot( geometry.corner(i) - yk );
                    if (under < 0.0){
                        continue;
                    }
                    alpha = std::min(alpha, (geometry.corner(i) - iterate).two_norm2() - maxDist*maxDist)
                    / ( 2.0 * under );
                }
            }
            if (alpha >= 1.0){
                return maxDist;
            }else {
                iterate = iterate + alpha*(yk - iterate);
            }

            k+=1;
        }

        //std::cout << iterate << std::endl;
        return maxDist;
    }

    template<class GridView>
    static void checkMesh(const GridView& gridView){
        Tutorial::StatisticsAccumulator meshSizeStatistic{};
        Tutorial::StatisticsAccumulator meshQualityStatistic{};
        double numElements = (double) gridView.size(0);

        for (const auto& element : elements(gridView)){
          auto geom = element.geometry();

          // Calculate the diameter
          double t = Tutorial::diameter(geom);
          double v = Tutorial::inscribedRadius(geom);

          meshQualityStatistic.add_datum(v/t);
          meshSizeStatistic.add_datum(t);

          // std::cout << element.type() << " with Diameter " << t << " " << v << std::endl;
        }

        std::cout << "Mesh Statistics [Min, Max] with Mean +- Variance:" << std::endl 
                  << "-------------------------------------------------" << std::endl
            << "Mesh Shape Regularity = [" 
            << meshQualityStatistic.calculate_min() << ", " << meshQualityStatistic.calculate_max() << "] with "
            << meshQualityStatistic.calculate_mean() << " +- " << meshQualityStatistic.calculate_std_deviation() << std::endl
            << "Mesh Diameter         = [" 
            << meshSizeStatistic.calculate_min() << ", " << meshSizeStatistic.calculate_max() << "] with "
            << meshSizeStatistic.calculate_mean() << " +- " << meshSizeStatistic.calculate_std_deviation() << std::endl;

    }


}