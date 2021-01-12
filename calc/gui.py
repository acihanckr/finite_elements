from tkinter import ttk, Tk, Text, BOTH, Label, messagebox, VERTICAL, HORIZONTAL, BOTTOM
import sympy
import numpy as np
from sympy import pprint, exp, diff, latex, Symbol, collect, expand, preview
from sympy.printing.latex import LatexPrinter, print_latex
from sympy.core.function import UndefinedFunction, Function
from sympy.abc import x,y,epsilon
from PIL import Image, ImageTk
from io import StringIO, BytesIO
from IPython.lib.latextools import latex_to_png
import sys
from expansion import *


class MyLatexPrinter(LatexPrinter):
    """Print derivative of a function of symbols in a shorter form.
    """
    def _print_Derivative(self, expr):
        function, *vars = expr.args
        if not isinstance(type(function), UndefinedFunction) or \
           not all(isinstance(i, Symbol) for i, degree in vars):
            return super()._print_Derivative(expr)

        # If you want the printer to work correctly for nested
        # expressions then use self._print() instead of str() or latex().
        # See the example of nested modulo below in the custom printing
        # method section.
        return "{{{}}}_{{{}}}".format(
            self._print(Symbol(function.func.__name__)),
                        ''.join(self._print(i) for i, degree in vars for j  in range(degree)))

    def print_my_latex(self, expr):
        """ Most of the printers define their own wrappers for print().
        These wrappers usually take printer settings. Our printer does not have
        any settings.
        """
        print(self.doprint(expr))



class Layout:


    def __init__(self, master):

        #create the latex printer
        self.printer = MyLatexPrinter()
        #setup the main window, with height  = 500, width = 200
        self.master = master
        self.master.geometry('640x480')

        #setup ttk widget styles
        self.style = ttk.Style()
        self.style.configure('TFrame', background='wheat')
        self.style.configure('TLabel', background='wheat')

        self.split = ttk.PanedWindow(self.master, orient = HORIZONTAL)
        self.split.pack(fill = BOTH, expand = True)

        #create the main frame to hold all other widgets for layout
        self.main_frame = ttk.Frame(master)
        #self.main_frame.pack(fill = BOTH, expand = True)

        self.frame1 = ttk.Frame(self.split)
        self.frame2 = ttk.Frame(self.split)

        self.split.add(self.frame1)
        self.split.add(self.frame2)

        #create text input and position it to grid 0x0, set event bind to key release
        self.input = Text(self.frame1, width = 50, height = 10)
        self.input.config(wrap = 'word')
        self.input.pack(fill = BOTH, expand = True)

        #create text output, connect it to stdout,  and position it to grid 1x0 
        self.output = ttk.Label(self.frame1)
        self.output.pack(fill = BOTH, expand = True)

        #create calculation frame to hold button, functions and degree entries
        self.calc_frame = ttk.Frame(self.frame2)
        self.calc_frame.pack(fill = BOTH, expand = True)
        


        #create function entries, label and pack to calculation frame
        self.funcs_label = ttk.Label(self.calc_frame, text = 'Functions to Expand').pack()
        self.funcs_ent = ttk.Entry(self.calc_frame, width = 20)
        self.funcs_ent.pack()


        #create degree entry and place to grid 0x1
        self.degree_label = ttk.Label(self.calc_frame, text = 'Degree of Expansion').pack()
        self.degree_ent = ttk.Entry(self.calc_frame, width = 20)
        self.degree_ent.pack()
                
        #create result button and place to grid 1x1 on calculation frame
        self.result_button = ttk.Button(self.calc_frame, text = 'Solve', command = self.find_result)
        self.result_button.pack()

        #create the preview button and place to grid 1x1 on calculation frame
        self.preview_button = ttk.Button(self.calc_frame, text = 'Preview', command = self.preview_exp)
        self.preview_button.pack()


        #create result frame
        self.result_frame = ttk.Frame(self.frame2)
        self.result_frame.pack(fill = BOTH, expand = True, side = BOTTOM, anchor = 's')

        #create the result text output and position to 1x1
        self.result_text = Text(self.result_frame, width = 50, height = 15)
        self.result_text.pack(fill = BOTH, expand = True)
        #scroll bar for result text
        xscroll = ttk.Scrollbar(self.result_frame, orient = HORIZONTAL, 
                    command = self.result_text.xview)
        xscroll.pack(fill = x, expand = True, anchor = 'n')
        self.result_text.config(xscrollcommand = xscroll.set)
        
        # #set up column and row weights for resizing
        # self.main_frame.rowconfigure(0, weight = 1)
        # self.main_frame.rowconfigure(1, weight = 1)
        # self.main_frame.columnconfigure(0, weight = 1)
        # self.main_frame.columnconfigure(1, weight = 1)

    def make_image(self,exp):
        old_out = sys.stdout
        new_out = StringIO()
        sys.stdout = new_out
        self.printer.print_my_latex(exp)
        exp = '$' + new_out.getvalue().strip('\n') + '$'
        sys.stdout = old_out
        obj = BytesIO()
        preview(exp, viewer='BytesIO', output = 'png', dvioptions=['-D','300'], outputbuffer=obj)
        obj.seek(0)
        img = Image.open(obj).convert("RGBA")
        imgnp = np.array(img)
        white = np.sum(imgnp[:,:,:3], axis=-1)
        white_mask = np.where(white == 255*3, 1, 0)
        imgnp[:,:,-1] = np.where(white_mask, 0, imgnp[:,:,-1])
        img = ImageTk.PhotoImage(Image.fromarray(np.uint8(imgnp)))
        obj.close()
        return img

    def preview_exp(self):
        u1 = Function('u1')(x,y,z)
        u2 = Function('u2')(x,y,z)
        exp = self.input.get('1.0', 'end')
        exp = parse_expr(exp, locals())
        self.output.img = self.make_image(exp)
        self.output.config(image = self.output.img)


    def find_result(self):

        #get entries
        exp = self.input.get('1.0', 'end')
        funcs = self.funcs_ent.get().split(', ')
        degree = self.degree_ent.get()
        solver = Solver(exp, funcs, int(degree))
        self.result_img = self.make_image(solver.solve())
        self.result_text.image_create('1.0', image = self.result_img)

        
        

def main():
    root = Tk()
    root.title('Expansion Calculator')
    root.geometry('640x480')
    gui = Layout(root)
    root.mainloop()


if __name__ == '__main__': main()