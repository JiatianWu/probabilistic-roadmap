This project is developed by Jiatian Wu (jiatian1) and Jucheng Zhang (jucheng1). 

1. How to Run
First run " mex planner " in Matlab to compile. To run Delaunay Neighborhood Connection based Probabilistic Roadmap, run " main(1) ". To run original Nearest Neighborhood Connection based PRM, run " main(0)".

Current arm configuration is same as hw2. Default start is [pi/2 pi/4 pi/2 pi/4 pi/2], default goal is [pi/8 3*pi/4 pi 0.9*pi 1.5*pi], default map is map1.txt.

2. Parameters Setting
Default parameters including neighborhood radius size and sampling points number are set on planner.cpp. For Delaunay Neighborhood Connection, set radius on line 1021 and sampling points number on line 1022. For Nearest Neighborhood Connection, set radius on line 941 and sampling points number on line 942.

To choose optimal parameters, please look at radius_map1.txt to choose a sampling points number and corresponding radius. This lookup table is constructed beforehand. To see how to constructed the radius lookup table, please install the gurobi solver first(http://www.gurobi.com/downloads/gurobi-optimizer). Then, you should add the code find\_delaunay\_c++.cpp from \get\_lookup\_table to the installed gurobi folder, and modify the Makefile to compile it. Plus, you need to modify the 716th line in find\_delaunay\_c++.cpp for the location of input map, and the 719th line for the location of output lookup table. 

Note: It might be very difficult to install gurobi and complie the code, please feel free to send emails to Contact for any questions.

3. Visualization
In addition to the path generated from start to end as hw2, the output in the command windows will also includes several metrcis: edges generated, average neighbors per node, components of the genrated graph, time taken. 

4. Contact
If you have problems when compiling or run the code, please contact us: Jiatian Wu (jiatian1@andrew.cmu.edu) and Jucheng Zhang (jucheng1@andrew.cmu.edu).
