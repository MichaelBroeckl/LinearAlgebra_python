# Students-Oriented Linear Algebra Toolkit

### Youtube Video: https://youtu.be/MFPfyqXcsEU

## Description:

This project provides a students-friendly framework experimenting with matrix operations 
and 3D vector geometry using Python. It includes Gaussian elimination, determinant computation
and geometric applications of vector math such as cross products, angles and distances.
The goal is to help students gain a deeper understanding of Linear Algebra and how to combine
it with Python programming. Therefore the script is extensively documented. To ensure the functions are not changed in a destructive way, unit tests were added.

## Project Files

- `project.py`: The main Python file. It includes all the required classes, methods and
functions necessary for functionality as well as a main function for interactive demonstration 
of Linear Algebra calculations
- `test_project.py`: Contains unit tests for the three core functions using pytest.
- `README.md`: You're reading it. This file explains the purpose, structure, and usage of the project.

## Classes

#### Matrix Class

This class is used to create a matrix out of a list of lists, calculate it's row-echelon form 
and from that, the determinant.

 #### Vector Class

 This class is used to create 3D-Vectors and use the following essential Linear Algebra operations on them:
 - Dot Product
 - Cross product
 - Magnitude 
 - Angle between vectors
 - Projection onto another vector
 - Scaling by a scalar

## Functions

- `parallelepiped_volume(v1, v2, v3)`: Computes the volume spanned by three vectors using the gaussian elimination and the determinant.
- `distance_point_to_line(point, line_point, direction)`: Uses the cross product and the magnitude methods to find the shortest distance from a point to a line in 3D Space.
- `distance_point_to_plane(point, plane_point, v1, v2)`: Uses the cross product and the magnitude methods to find the shortest distance from a point to a plane in 3D Space.

## Demonstration:

The `main()` function shows how all features can be used and customized. Students can modify vectors and matrices to become familiar with the concepts and maybe even check their homework.

## Design Considerations:

Due to flexibility, clarity and structure I used object oriented programming.
For the project, three functions were required and I wasn't sure if methods would be counted as functions by the automated CS50 submission test. Otherwise I maybe would have stuck to OOP except for the `main()` function obviously.

