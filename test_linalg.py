from math import isclose
from linalg import parallelepiped_volume, distance_point_to_line, distance_point_to_plane

def test_parallelepiped_volume():
    v1 = [1, 0, 0]
    v2 = [0, 1, 0]
    v3 = [0, 0, 1]
    assert parallelepiped_volume(v1, v2, v3) ==1

def test_distance_point_to_line():
    point = (1, 2, 3)
    line_point = [1, 0, 0]
    direction = [0, 1, 0]
    result = distance_point_to_line(point, line_point, direction)
    assert isclose(result, 3.0, abs_tol=1e-5)

def test_distance_point_to_plane():
    point1 = [0, 0, 1]
    plane_point = [0, 0, 0]
    vector1 = [1, 0, 0]
    vector2 = [0, 1, 0]
    result = distance_point_to_plane(point1, plane_point, vector1, vector2)
    assert isclose(result, 1.0, abs_tol=1e-5)

def main():
    test_parallelepiped_volume()
    test_distance_point_to_line()
    test_distance_point_to_plane()

if __name__ == '__main__':
    main()