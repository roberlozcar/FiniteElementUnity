# Finite Element Solver for Unity


## Introduction

We implemented a tridimensional finite elements solver for an arbitrary object and the collisions with planes and spheres.
		
The collision objects are considered planes unless you use the tag "Sphere" in the unity editor. Also the spheres must have a sphere collider with the correct radius.
		
We need the geometry of the object, in .obj format, and the volume discretization in tetrahedrons (we suggest using TetGen to calculate it), the list of nodes in .node format and the tetrahedrons list in .ele. These file must be in the folder Assets/Resources. Furthermore, in the unity editor, you have to change the "Ident" field of the object to match the name of these files. There are two examples in the Resources folder, an sphere and the Stanford bunny.
	
To start the simulation, you have to unmark the pause field in physic manager object editor.


## Options

You can modify the mechanical parameters of the object and the gravity acceleration in the unity editor of the object.
		
You can choose the solver of the physical system: Symplectic, which is faster but it need lower steps and it can be unstable in the collision and implicit, which is slower but it can use higher step and it is stable in collision (if the penetration is not big). This method can use dense or sparse matrix, we recommend sparse because the solver is faster.
	
Also you can choose what geometry, the tetrahedron nodes or the mesh vertex, is used to compute the collision. The mesh vertex produce a more realistic contact but it is slower.
		
You can also add more collision objects and load more objects to simulate, with a computational cost, adding more elements to the coll obstacles and loaded objects field in physic manager object.