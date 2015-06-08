#!/bin/sh


#python gprof2dot.py -f pstats output.pstats | dot -Tpng -o output.png



python -m cProfile -o output.pstats main_example.py
#python gprof2dot.py -n0.1 -e0.1 -z "scripts:95:stoner0K_rotation" -f pstats output.pstats | dot -Tpng -o output.png

#python gprof2dot.py -n0.1 -e0.1 -f pstats output.pstats | dot -Tpng -o output.png
