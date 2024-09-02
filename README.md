# Tube2FEM - 
Tube2FEM is a general-purpose highly-automated pipeline for flow related processes in (embedded) tubular objects

## Getting started

[Click here](https://www.gibboncode.org) 

![Tetrahedral mesh](https://www.gibboncode.org/img/bunnyMesh.gif) 

The image doesn't work


Authors
-------
Tube2FEM is developed by:

  * Hani Cheikh Sleiman
  * Shiyu Wang

Licence
-------
turtleFSI is licensed under the GNU GPL, version 3 or (at your option) any
later version. turtleFSI is Copyright (2016-2019) by the authors.

Documentation
-------------
For an introduction to Tube2FEM, and tutorials, please refer to the [documentation](https://readthedocs.org/).

If you wish to use Tube2FEM for journal publications, please refer to the [JOSS publication](https://joss.theoj.org/papers/10.21105/joss.02089#):

Cheikh Sleiman et al., (2024). Tube2FEM: Tube2FEM: a general-purpose highly-automated pipeline for flow related processes in (embedded) tubular objects. Journal of ... https://...


Installation
------------



Use
---
Run turtleFSI with all the default parameters::
   ``turtleFSI``

See all the command line parameters run the following command::
  ``turtleFSI -h``

Run a specific problem file::
  ``turtleFSI --problem [path_to_problem]``

When calling a specific problem file, turtleFSI will first look for the file name locally, then check if the file name is present in the directory "/turtleFSI/problems/".
Please refere to the [documentation](https://turtlefsi2.readthedocs.io/en/latest/) to learn how to define a new problem file and for a more complete description of usage.


Contact
-------
The latest version of this software can be obtained from

  https://github.com/KVSlab/turtleFSI

Please report bugs and other issues through the issue tracker at:

  https://github.com/KVSlab/turtleFSI/issues
