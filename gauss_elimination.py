class Matrix:
    """ 
    Expects the dimensions of the coefficient matrix m as the number of rows and n as the number of columns as int and
    A, which can be an coefficient- or extended coefficient matrix in the form of list of lists and initializes
    the indices k=0, i=0 and j=0 and the variable for counting the rowswaps l=0.

    """

    def __init__(self, m, n, A):
        self.step = 0
        self.l = 0
        self.k = 0
        self.i = 0
        self.j = 0
        self.m = m
        self.n = n
        self.A = A
        


    def gauss_elimination(self):
        self.A = [[float(value) for value in row] for row in self.A]
        while True:
            self.find_smallest_k()
            if self.k == None:
                self.j +=1 
                if self.j > self.n:
                    return self.A
                else:
                    continue
            if self.k != self.i:
                self.A[self.k], self.A[self.i] = self.A[self.i], self.A[self.k]
                self.l += 1
            self.row_eliminations()
            self.i += 1
            self.j += 1
            self.step += 1
            print(f'Step {self.step} :')
            for row in self.A:
                print(row)
            if self.i >= self.m - 1 or self.j > self.n - 1:
                return self.A, self.l

    def find_smallest_k(self):
        self.k = None
        for k in range(self.i, self.m):
            if self.A[k][self.j] != 0:
                self.k = k
                break
        
    def row_eliminations(self):
        for k in range(self.i + 1, self.m):
            if self.A[k][self.j] != 0:
                factor = self.A[k][self.j] / self.A[self.i][self.j]
                for x in range(self.j, self.n):
                    self.A[k][x] += self.A[self.i][x] * (-1) * factor 

    

    

def main():
    A = [
            [1, 0, 1, -3],
            [2, 0, 2, 1],
            [0, -2, -1, 0],
            [-1, 1, 3, 1],
        ]
    MatrixA = Matrix(4,4,A)
    Ã = MatrixA.gauss_elimination()
    print(f'With l = {Ã[1]} rowswaps, \nthe final matrix in echelon form is:')
    for row in Ã[0]:
        print(row)



if __name__ == "__main__":
    main()