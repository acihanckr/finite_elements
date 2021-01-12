import sympy
from sympy import *
from sympy.printing.latex import LatexPrinter, print_latex
from sympy.core.function import Function
from sympy.parsing.latex import parse_latex
from sympy.parsing.sympy_parser import parse_expr
from sympy.abc import x,y,z,epsilon
from sys import stdout
from gui import *

class Solver:


    def __init__(self, exp, funcs, degree):
        self.exp = exp
        self.fs = funcs
        self.degree = degree
        self.rep = {f:'self.lst[{}]'.format(i) for f,i in zip(self.fs,range(len(self.fs)))}
        self.lst = []
        for u in funcs:
            us = []
            for i in range(degree):
                name = u + str(i)
                f = Function(name)(x,y,z)
                us.append(f)
            genu = us[0]
            for d, uin in enumerate(us[1:]):
                genu = genu + epsilon**(d+1)*uin
            self.lst.append(genu)    
        self.my_parse_exp()

        
    def my_parse_exp(self):
        self.exp = self.exp.replace(' ', '')
        for key, value in self.rep.items():
            self.exp = self.exp.replace(key, value)
    
    def solve(self):
        self.result = eval(self.exp)
        self.result = collect(expand(self.result),epsilon, evaluate=False)
        return self.result


def main():
    exp = 'diff(u1, x, x)+ diff(u1, y, y)/epsilon**2+ diff(u1,x,x)/epsilon + diff(u2, z)/epsilon**2 - 2*(u1**2+u2**2-1)*u1/epsilon**2'
    solver = Solver(exp, ['u1', 'u2'], 3)
    pprint(solver.solve())


if __name__ == '__main__': main()