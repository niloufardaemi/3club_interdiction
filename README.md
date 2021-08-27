# Code for interdicting the maximum 2-club in a graph

This code accompanies the paper "Interdicting Low-Diameter Cohesive Subgroups in
Large-Scale Social Networks" and is written in C++. If you wish to use or cite this code, please cite the paper:


@article{Niloufar2020InterClubs,
	author = {Niloufar Daemi and Juan S. Borrero and Balabhaskar Balasundaram},
	journal = {INFORMS Journal on Optimization},
	month = {August},
	note = {Accepted for publication.},
	title = {Interdicting low-diameter cohesive subgroups in large-scale social networks},
	year = {2021}
	}
  
  
This repository includes two folders:

1. exact_separation: used to solve the separation problem (maximum s-club) to optimality.
2. heuristic_separation: used to first try heuristic approaches to solve the separation problem and solving it to optimality only if needed.


# Compiling the code

The following steps show how to compile and run the code to find the maximum 2-club in the interdicted graph in a Windows environment.

1. Download or clone the repository to your machine.
2. Build the project.
3. In the command prompt, change the directory to the one that the project is stored and type the following command with space instead of each comma: project name.exe, dataset format (dimacs/snap_d), dataset directory/, dataset name, value of alpha, and hit enter.



# Terms and Use:

MIT License

Copyright (c) 2021 Niloufar Daemi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



# Acknowledgments

We would like to thank Austin Buchanan, Hosseinali Salemi and Hamidreza Validi for providing
us the codes used in [Salemi and Buchanan 2020](https://link.springer.com/article/10.1007/s12532-020-00175-6) and [Validi and Buchanan 2020](https://pubsonline.informs.org/doi/abs/10.1287/ijoc.2019.0914) and freely offering their help
to integrate it with our code. All the files in our repository except the main.cpp are part of therir codes with some modifications. Their codes in the original form respectively availble at [1](https://github.com/halisalemi/ParsimoniousKClub) and [2](https://github.com/hamidrezavalidi/The-Optimal-Design-of-Low-Latency-Virtual-Backbones).
