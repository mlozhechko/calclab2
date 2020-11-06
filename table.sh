#!/bin/bash

echo ""
echo "fpi (T = 0.008)"
echo ""
echo "ordinary"
echo ""
./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method fpi -criteria ordinary -solution ./data/sys2/12/solution.dat -eps 1E-7
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method fpi -criteria ordinary -solution ./data/sys2/12/solution.dat -eps 1E-4

echo ""
echo "solution"
echo ""
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method fpi -criteria solution -solution ./data/sys2/12/solution.dat -eps 1E-7
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method fpi -criteria solution -solution ./data/sys2/12/solution.dat -eps 1E-4

echo ""
echo "delta"
echo ""
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method fpi -criteria delta -solution ./data/sys2/12/solution.dat -eps 1E-7
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method fpi -criteria delta -solution ./data/sys2/12/solution.dat -eps 1E-4

echo "jacobi"
echo ""
echo "ordinary"
echo ""
./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method jacobi -criteria ordinary -solution ./data/sys2/12/solution.dat -eps 1E-7
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method jacobi -criteria ordinary -solution ./data/sys2/12/solution.dat -eps 1E-4

echo ""
echo "solution"
echo ""
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method jacobi -criteria solution -solution ./data/sys2/12/solution.dat -eps 1E-7
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method jacobi -criteria solution -solution ./data/sys2/12/solution.dat -eps 1E-4

echo ""
echo "delta"
echo ""
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method jacobi -criteria delta -solution ./data/sys2/12/solution.dat -eps 1E-7
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method jacobi -criteria delta -solution ./data/sys2/12/solution.dat -eps 1E-4

echo ""
echo "seidel"
echo ""
echo "ordinary"
echo ""
./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method seidel -criteria ordinary -solution ./data/sys2/12/solution.dat -eps 1E-7
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method seidel -criteria ordinary -solution ./data/sys2/12/solution.dat -eps 1E-4

echo ""
echo "solution"
echo ""
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method seidel -criteria solution -solution ./data/sys2/12/solution.dat -eps 1E-7
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method seidel -criteria solution -solution ./data/sys2/12/solution.dat -eps 1E-4

echo ""
echo "delta"
echo ""
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method seidel -criteria delta -solution ./data/sys2/12/solution.dat -eps 1E-7
 ./calc2 -precision double -matrix ./data/sys2/12/matrix.dat -vector ./data/sys2/12/vector.dat\
 -method seidel -criteria delta -solution ./data/sys2/12/solution.dat -eps 1E-4