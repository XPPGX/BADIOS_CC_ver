/*
---------------------------------------------------------------------
 This file is part of BADIOS framework
 Copyright (c) 2012,
 By:    Ahmet Erdem Sariyuce,
        Erik Saule,
        Kamer Kaya,
        Umit V. Catalyurek
---------------------------------------------------------------------
 For license info, please see the README.txt and LICENSE.txt files in
 the main directory.
---------------------------------------------------------------------
*/

======
BADIOS
======

- To compile the code, make sure you have a version of gcc and type
    $> MODE=gcc OPT=yes PRINT_BC_VALUES=yes make -e

- Extra compilation options are available in makefile.in with
  explanations

- To run give the following arguments:
    <input filename> 
    <nTry>
    <0 for BADIOS, 1 for base-brandes> // Ordering is enabled 
                                          as default in BADIOS
    <side-vertex_disable>              // 1 to disable,
                                          0 to enable side vertex removal
    <articulation-vertex_disable>      // 1 to disable,
                                          0 to enable articulation vertex copy
    <bridge_disable>                   // 1 to disable,
                                          0 to enable bridge removal
    <degree1-vertex_disable>           // 1 to disable,
                                          0 to enable degree-1 vertex removal
    <identical-vertex_disable>         // 1 to disable,
                                          0 to enable identical vertex removal 
                                          (both type-I and type-II)

- For <input filename>, Chaco and Matrix Market formats are supported. Chaco
  format input files should have an extension of ".graph".
  Unsymmetic Matrix Market files can be given as input since there is a
  symmetrization procedure for them. However, Chaco files should
  be given as symmetric.

- Example run where articulation vertex copy and degree-1 vertex techniques
  are enabled and others are disabled:
     $> ./bc-seq-brandes power.graph 1 0 1 0 1 0 1"

- For any question or problem, please contact:
    aerdem@bmi.osu.edu
    esaule@bmi.osu.edu
    kamer@bmi.osu.edu
    umit@bmi.osu.edu
