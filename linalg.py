"""
Author: Michael Broeckl, GitHub: MichaelBroeckl edX: MB_2501_43J5
Course: CS50P, CS50's Introduction to Programming with Python
Date: April 22, 2025
City: Salzburg 
Country: Austria

Students-Oriented Linear Algebra Toolkit

This script provides a students-friendly framework experimenting with matrix operations 
and 3D vector geometry using Python. It includes Gaussian elimination, determinant computation
and geometric applications of vector math such as cross products, angles and distances.
The goal is to help students gain a deeper understanding of Linear Algebra and how to combine
it with Python programming. Therefore the script is extensively documented.
"""

from math import sqrt, acos, degrees


# -------------------- MATRIX CLASS ------------------

class Matrix:
    """ 
    A class to perform Gaussian elimination and compute the determinant of a matrix.

    This class supports operations on both standard matrices (A) and extended 
    matrices (A|b) (augmented with a column vector, e.g., for systems of equations).

    Attributes:
        A (list of lists of floats): The matrix.
        m (int): Number of rows.
        n (int): Number of columns (excluding the extended column if applicable)
        extended (bool): Indicates whether the matrix is extended (i.e., includes a results_column)
        i, j, k (int): Indices used in elimination algorithm.
        l (int): Count of row swaps (used for determinant sign).
        step (int): Step counter for display purposes.
        det_sum (float): Running product of diagonal elements for determinant calculation.

    Args:
    A(list of lists of numbers): Coefficient matrix (can include results_column).
    extended (bool): Whether A is an extended matrix (A|b).
    show_calculation (bool, optional): Whether the whole calculation is printed (defaults to False).
    """

    def __init__(self, A, extended: bool, show_calculation=False):
        self.det_sum = 1
        self.step = 1
        self.l = 0
        self.k = 0
        self.i = 0
        self.j = 0
        self.A = A
        self.m = len(self.A)
        self.extended = extended
        self.show_calculation = show_calculation
        if self.extended:
            self.n = len(self.A[0]) - 1
        else:
            self.n = len(self.A[0])



    def gauss_elimination(self):
        """
        Performs Gaussian elimination to transform the matrix into row echelon form.
        
        Returns:
            tuple: (Echelon matrix as list of lists, number of row swaps, determinant)
        """
        self.A = [[float(value) for value in row] for row in self.A]
        while True:
            self.find_smallest_k()
            if self.k is None:
                self.j +=1
                if self.j > self.n:
                    return self.A, self.l, self.det()
                else:
                    continue
            if self.show_calculation:
                print(f'Step {self.step} :')
            self.row_swap()
            self.row_eliminations()
            self.i += 1
            self.j += 1
            self.step += 1
            if self.show_calculation:
                for row in self.A:
                    print(row)
            if self.i >= self.m - 1 or self.j > self.n - 1:
                return self.A, self.l, self.det()

    def find_smallest_k(self):
        """
        Finds the number of the first row >= i in the current column j where the value 
        of the element != 0.
        """
        self.k = None
        for k in range(self.i, self.m):
            if self.A[k][self.j] != 0:
                self.k = k
                break

    def row_eliminations(self):
        """
        Eliminates the element of the k-th row in the j-th column by adding the i-th row 
        multiplied by (-1) * (pivot-k-th row / pivot-i-th row) to the k-th row.
        """
        for k in range(self.i + 1, self.m):
            if self.A[k][self.j] != 0:
                factor = (-1) * self.A[k][self.j] / self.A[self.i][self.j]
                if self.show_calculation:
                    print(f'Row{k}new = Row{k} + Row{self.i} x {factor}')
                if self.extended:
                    for x in range(self.j, self.n + 1):
                        self.A[k][x] += self.A[self.i][x] * factor
                else:
                    for x in range(self.j, self.n):
                        self.A[k][x] += self.A[self.i][x] * factor

    def det(self):
        """
        Computes the determinant from the row echelon form.
        
        Returns:
            float: Determinant value.
        """
        number_of_nonzero_rows = sum(not all(x == 0 for x in sublist[:self.n])
                                     for sublist in self.A)
        for i in range(number_of_nonzero_rows):
            self.det_sum *= self.A[i][i]
        return (-1) ** self.l * self.det_sum

    def row_swap(self):
        """
        Swaps the respective rows if, in the j-th column, k != i.
        """
        if self.k != self.i:
            if self.show_calculation:
                print(f'swap Row{self.k} with Row{self.i}')
            self.A[self.k], self.A[self.i] = self.A[self.i], self.A[self.k]
            self.l += 1

# -------------------- VECTOR CLASS ------------------

class Vector:
    """ 
    A class for basic 3D vector operations including dot product, cross product,
    magnitude, angle and projection.
    
    Args:
        values (list [float]): A 3-dimensional vector.
    """

    def __init__(self, values):
        if len(values) != 3:
            raise ValueError('Only 3D vectors are supported.')
        self.values = values

    def __str__(self):
        return f'({', '.join(f"{v:.2f}" for v in self.values)})T'

    def dot(self, other):
        """ Computes the dot product between this vector and another.
        
        Args:
            other (Vector): The other vector.
            
        Returns:
            float: Scalar result of the dot product.
        """
        return sum(a*b for a, b in zip(self.values, other.values))

    def cross(self, other):
        """
        Computes the cross product between this vector and another.
        
        Args:
            other (Vector): The other vector.
        
        Returns:
            Vector: A new vector orthogonal to both."""
        a, b = self.values, other.values
        return Vector([
            a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0],
        ])

    def magnitude(self):
        """
        Computes the Euclidean norm (length) of the vector.
        
        Returns:
            float: Magnitude of the vector.
        """
        return sqrt(sum([a**2 for a in self.values]))

    def angle_with(self, other):
        """
        Computes the angle (in degrees) between this vector and another.
        
        Args:
            other (Vector): The other vector.
        
        Returns:
            float: Angle in degrees.
            
        Raises:
            ValueError: If either vector is zero-length.
        """
        magnitudes = self.magnitude()*other.magnitude()
        if magnitudes == 0:
            raise ValueError('Vector with magnitude 0 not allowed.')
        cos_theta = max(min(self.dot(other)/magnitudes, 1), -1)
        return degrees(acos(cos_theta))

    def scale(self, scalar):
        """
        Scales the vector by a given scalar.
        
        Args:
            scalar(float): Scalar multiplier.
            
        Returns:
            Vector: Scaled vector.
        """
        return Vector([scalar * a for a in self.values])

    def projection_onto(self, other):
        """
        Projects this vector onto another vector using the normalized direction of the other.
        
        Formula:
            proj_b(a) = (a . b ) * bn
            where bn = b / ||b|| is the normalized vector of b
        
        Args:
            other (Vector): The vector to project onto.
            
        Returns:
            Vector: The projected vector of self onto other.
        """

        unit = other.scale(1 / other.magnitude())
        return unit.scale(self.dot(unit))

# ---------------------- FUNCTIONS ------------------
def allow_vectors(*args):
    """Allows you to use vector objects or lists as vectors."""
    return [v.values if isinstance(v, Vector) else v for v in args]

def parallelepiped_volume(v1, v2, v3):
    """
    Calculates the volume of a parallelepiped defined by three 3d vectors
    using the absolute value of the determinant of the matrix composed of the vectors as columns.

    Args:
        v1, v2, v3 (list[float] or Vector)
        
    Returns:
        float: The volume of the parallelepiped.
    """
    v1, v2, v3 = allow_vectors(v1, v2, v3)
    matrix_data = [[v1[i], v2[i], v3[i]] for i in range(3)]
    matrix = Matrix(matrix_data, extended=False, show_calculation=False)
    det = matrix.gauss_elimination()[2]
    if det == 0:
        raise ValueError('The given vectors are not linearly independent'
                        'and therefore don\'t span a parallelepiped')
    return abs(det)

def distance_point_to_line(point, line_point, direction):
    """
    Calculates the shortest distance from a point to a line in 3d.
    
    The line is defined by a point on the line and a direction vector.
    This uses the formula (x means the cross product):
        distance = ||(point - line_point ) x direction|| / ||direction||
    
    Args:
        point (list[float] or Vector): Coordinates of the external point [x, y, z]
        line_point(list[float] or Vector): A point on the line [x, y, z]
        direction (list[float] or Vector): Direction vector of the line [x, y, z]

    Returns:
        float: The shortest distance from the point to the line as scalar.
    """
    point, line_point, direction = allow_vectors(point, line_point, direction)
    p = Vector(point)
    q = Vector(line_point)
    v = Vector(direction)
    pq = Vector([a-b for a, b in zip(p.values, q.values)])
    cross = pq.cross(v)
    return cross.magnitude() / v.magnitude()

def distance_point_to_plane(point, point_plane, vector1, vector2):
    """
    Calculates the shortest distance from a point to a plane in 3d.
    The plane is given by a point inside the plane and two vectors.
    The distance is calculated by the Formula:
    distance = |dot product of pq and n| / ||n||
    where:
        pq = vector point to point_plane
        n = norm vector of the plane
    Args:
        point(list[float] or Vector): Coordinates of a point in 3d space [x, y, z]
        point_plane(list[float] or Vector): Coordinates of a point inside the plane [x, y, z]
        vector1, vector2 (list[float] or Vector): Vectors that span the plane [x, y, z]
    Returns:
        float: The shortest distance from the point to the plane as scalar.
    """
    point, point_plane, vector1, vector2 = allow_vectors(point, point_plane, vector1, vector2)
    p = Vector(point)
    q = Vector(point_plane)
    pq = Vector([a-b for a,b in zip(p.values, q.values)])
    v1 = Vector(vector1)
    v2 = Vector(vector2)
    n = v1.cross(v2)
    return abs(pq.dot(n))/n.magnitude()



# ------------------------ MAIN -------------------------



def main():
    """
    Run demonstration of matrix and vector operations for educational/research use.
    Modify values below to explore different configurations.
    """

    # ------- Gauss Elimination Application -------
    # Edit Matrix A:
    A = [
            [1, 0, 1, -3],
            [2, 0, 2, 1],
            [0, -2, -1, 0],
            [-1, 1, 3, 1],
        ]

    # Add the argument "show_calculation=True" to see each step of the calculation.
    matrix = Matrix(A, extended=False, show_calculation=True)
    echelon_matrix, swaps, determinant = matrix.gauss_elimination()
    print('This is the final matrix in echelon form:')
    for row in echelon_matrix:
        print(row)
    print(f'with l = {swaps} row swaps and\n'
          f'the determinant = {determinant}')

    # ------ Dot and Cross Product, Magnitude, Angle, Scaling and Projection ------
    # Edit Vectors v and w and the scalar by which you want to scale n:
    v = Vector([1,2,3])
    w = Vector([3,2,1])
    scalar = 5
    n = v.cross(w)
    print(f'The dot product v.w = {v.dot(w)}')
    print(f'The cross product n = v x w = {n}')
    print(f'The angle between v and w, phi = {v.angle_with(w):.2f} degree')
    print(f'The scaled vector ñ = {n.scale(scalar)}')
    print(f'v projected onto w gives us vp = {v.projection_onto(w)}')

    # ------ Volume Parallelepiped ------
    # Edit the vectors a, b and c
    a = Vector([1, 1, 1])
    b = Vector([2, 3, 3])
    c = Vector([1, -1, 2])
    print(f'The Volume of the parallelepiped Vp = {parallelepiped_volume(a, b, c)}')

    # ------ Distance Point to Line ------
    # Edit the point and the line_point and direction.
    point = [1, 3, 5]
    line_point = [3, 1, 2]
    direction = [1, -1, 7]
    print('The distance point to line dl = '
         f'{distance_point_to_line(point, line_point, direction):.2f}')

    # ------ Distance Point to Plane ------
    # Edit point1, point_plane, vector1 and vector2
    point1 = [1, 1, 1]
    point_plane = [-1, -1, -1]
    vector1 = [2, 0, 1]
    vector2 = [1, -2, 1]
    print('The distance point to plane dp ='
         f' {distance_point_to_plane(point1, point_plane, vector1, vector2):.2f}')


if __name__ == "__main__":
    main()
