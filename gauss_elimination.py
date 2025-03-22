class Matrix:
    """ 
    Expects the dimensions of the coefficient matrix: m as the number of rows 
    , n as the number of columns as int, A in the form of list of lists andinitializes
    the indices k=0, i=0 and j=0, the variable for counting the rowswaps l=0
    and the values for step=1 and det_sum=1.
    In case the Matrix is not an extended coefficient matrix set extended from the 
    default extended=True to extended=False.
    """

    def __init__(self, A, extended=True):
        self.det_sum = 1
        self.step = 1
        self.l = 0
        self.k = 0
        self.i = 0
        self.j = 0
        self.A = A
        self.m = len(self.A)
        self.extended = extended
        if self.extended:
            self.n = len(self.A[0]) - 1
        else:
            self.n = len(self.A[0])



    def gauss_elimination(self):
        self.A = [[float(value) for value in row] for row in self.A]
        if self.extended:
            print('(A|b)=')
        else:
            print('(A)=')
        for row in self.A:
                print(row)
        while True:
            self.find_smallest_k()
            if self.k == None:
                self.j +=1 
                if self.j > self.n:
                    return self.A, self.l, self.det()
                else:
                    continue
            print(f'Step {self.step} :')
            self.row_swap()
            self.row_eliminations()
            self.i += 1
            self.j += 1
            self.step += 1
            for row in self.A:
                print(row)
            if self.i >= self.m - 1 or self.j > self.n - 1:
                return self.A, self.l, self.det()

    def find_smallest_k(self):
        self.k = None
        for k in range(self.i, self.m):
            if self.A[k][self.j] != 0:
                self.k = k
                break
        
    def row_eliminations(self):
        for k in range(self.i + 1, self.m):
            if self.A[k][self.j] != 0:
                factor = (-1) * self.A[k][self.j] / self.A[self.i][self.j]
                print(f'Row{k}new = Row{k} + Row{self.i} x {factor}')
                if self.extended:
                    for x in range(self.j, self.n + 1):
                        self.A[k][x] += self.A[self.i][x] * factor
                else:
                    for x in range(self.j, self.n):
                        self.A[k][x] += self.A[self.i][x] * factor 
    
    def det(self):
        number_of_nonzero_rows = sum(not all(x == 0 for x in sublist[:self.n]) for sublist in self.A)
        for i in range(number_of_nonzero_rows):
            self.det_sum *= self.A[i][i]
        return (-1) ** self.l * self.det_sum
        
    def row_swap(self):
        if self.k != self.i:
                print(f'swap Row{self.k} with Row{self.i}')
                self.A[self.k], self.A[self.i] = self.A[self.i], self.A[self.k]
                self.l += 1

    

def main():
    A = [
            [1, 0, 1, -3],
            [2, 0, 2, 1],
            [0, -2, -1, 0],
            [-1, 1, 3, 1],
        ]
    B = [
        [4, -1],
        [3, -2],
    ]

    C = [
        [1, -2, 1, 1],
        [0, 1, 1, 1],
        [1, 0, 1, 0],
    ]
    MatrixA = Matrix(A, extended=False)
    Ã = MatrixA.gauss_elimination()
    print(f'This is the final matrix in echelon form,\n'
          f'with l = {Ã[1]} rowswaps and\n'
          f'the Determinant det(Ã) = {Ã[2]}')



if __name__ == "__main__":
    main()