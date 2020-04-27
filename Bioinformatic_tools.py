import tkinter
from tkinter import *
from collections import deque


def fasta_option():
    Label(window, text="Give there a FASTA file", font="none 12 bold").pack()
    entry = Entry(window, width=50)
    entry.pack()  
    Button(window, text="OK", command =lambda: fasta_worker(entry.get())).pack()
def fasta_subsequence(name_seq, second_dict, first, second):
    sequence = second_dict.get(name_seq)
    seq = sequence.splitlines()  
    sequence = ''.join(seq)  
    fir = int(first) - 1
    sec = int(second) - 1
    Label(window, text=sequence[fir:sec]).pack()
def fasta_output(num, first_dict, second_dict, sequences, value):
    clearFrame()
    for i in range (0,num):
        if (i == (value-1)):
            Label (window, text=sequences[i], font="none 12 bold").pack()
            Label (window, text="Description of the sequence", font="Helvetica 10 bold").pack()
            Label (window, text= first_dict.get(sequences[i])).pack()
            Label (window, text="Sequence:", font="Helvetica 10 bold").pack()
            Label (window, text=second_dict.get(sequences[i])).pack()
            Label (window, text="Length of sequence:", font="Helvetica 10 bold").pack()
            Label (window, text= (len(second_dict.get(sequences[i]))+1)).pack()
            Label (window, text= "Subsequence", font="Helvetica 10 bold").pack()
            Label (window, text= "Give two numbers between " + '0 and ' + str(len(second_dict.get(sequences[i]))+1), font="Helvetica 10 bold").pack()
            name = sequences[i]
            entry1 = Entry(window, width=30)
            entry2 = Entry(window, width=30)
            entry1.pack()
            entry2.pack()
            Button(window, text="OK", command =lambda: fasta_subsequence(name, second_dict, entry1.get(), entry2.get())).pack()      


def fasta_worker(e):
    sequences = []
    r = IntVar()
    num = 0
    first_dict = {} #description of sequence
    second_dict = {} #sequence
    with open (e) as f:
        for line in f:            
            if len(line.strip()) != 0:
                splitter = line.split()
                if (splitter[0][:1] == '>'):
                    num += 1
                    sequences.append(splitter[0][1:])
                    string = ' '.join(splitter[1:])                    
                    first_dict[sequences[num-1]] = string
                    second_dict[sequences[num-1]] = None
                else:
                    string = '\n'.join(splitter)
                    if (second_dict[sequences[num-1]] == None):
                        second_dict[sequences[num-1]] = string
                    else:
                        second_dict[sequences[num-1]] = second_dict[sequences[num-1]] + '\n' + string                    

    Label(window, text="Choose one of these sequences", font="none 12 bold").pack()
    num = 1            
    for sequence in sequences:
        Radiobutton(window, text=sequence, variable=r, value=num).pack()
        num += 1
    Button(window, text="OK", command =lambda: fasta_output(num,first_dict, second_dict, sequences, r.get())).pack()


def hamming_worker(first, second):
    res = 0
    if (len(first) == len(second)):
        for i in range (0, len(first)):
            if (first.upper()[i] != second.upper()[i]):
                res += 1
        Label(window, text="Hamming distance between these sequences is: " + str(res), font="none 12 bold").pack()
    else:
        Label(window, text="The two sequences don't have a same length", font="none 12 bold").pack()
    Button(window,text="BACK", command =lambda: hamming_option()).pack()

def hamming_choice(value):
    if (value == 1):
        fasta_option()

    else:
         Label(window, text="Write the first sequence:",font="none 12 bold").pack()
         entry1 = Entry(window, width=150); entry1.pack()         
         Label(window, text="Write the seqond sequence:",font="none 12 bold").pack()
         entry2 = Entry(window, width=150); entry2.pack()         
         Button(window, text="OK", command =lambda: hamming_worker(entry1.get(), entry2.get())).pack()

def hamming_option():
    clearFrame()
    Label (window, text="Choose if you want to load a FASTA file or only write two sequences:",font="none 12 bold").pack()
    r = IntVar()
    r1 = Radiobutton(window, text="FASTA file", variable=r, value=1).pack()
    r2 = Radiobutton(window, text="Write two sequences", variable=r, value=2).pack()
    Button(window, text="OK", command =lambda: hamming_choice(r.get())).pack()

def edit_backtracking(matrix, first, second):
    fir = []; sec = []
    lf = len(first); ls = len(second)    
    while (matrix[lf][ls] != 0):
        act = matrix[lf][ls]
        left = matrix[lf][ls-1]
        up = matrix[lf-1][ls]
        if (act == matrix[lf-1][ls-1] or (act-1) == matrix[lf-1][ls-1]):
            fir.append(first[lf-1]); sec.append(second[ls-1])
            lf -= 1; ls -= 1
        elif (act-1 == matrix[lf-1][ls]):
            fir.append(first[lf-1]); sec.append('-')
            lf -= 1
        elif (act-1 == matrix[lf][ls-1]):
            fir.append('-'); sec.append(second[ls-1])
            ls -= 1   
        elif (matrix[lf][ls] == 0): break
    if (lf == 2 and ls != 2):       
        nm = ls - lf
        for i in range (nm, 0, -1):            
            fir.append('-'); sec.append(second[ls-1])
            ls -= 1
    elif (lf != 2 and ls == 2):
        nm = lf - ls
        for j in range (nm, 0, -1):
            fir.append(first[lf-1]); sec.append('-')
            lf -= 1
    print(fir); print(sec)
    fir.reverse(); print(fir)
    first_string = ' '.join(fir)
    sec.reverse(); print(sec)
    second_string = ' '.join(sec)
    Label(window, text="Sequence alignment: ", font="none 12 bold").pack()
    Label(window, text=first_string).pack()
    Label(window, text=second_string).pack()    
    
def make_matrix(first, second):
    matrix = []
    for i in range(first):
        matrix.append([0]*second)
    return matrix

def edit_worker(first, second):
    matrix = make_matrix(len(first)+1, len(second)+1)
    num = 0
    #first row
    for i in range(len(second)+1):
        matrix[0][i] = num
        num += 1
    num = 0
    #first column
    for j in range(len(first)+1):
        matrix[j][0] = num
        num += 1
    for m in range (1,len(first)+1):
        for n in range (1, len(second)+1):
            if (first.upper()[m-1] == second.upper()[n-1]):
                matrix[m][n] = matrix[m-1][n-1]
            else:
                matrix[m][n] = min(matrix[m-1][n-1]+1, matrix[m][n-1]+1, matrix[m-1][n]+1)
        edit_distance = matrix[m][n]
    Label(window, text="Edit distance: " + str(edit_distance),font="none 12 bold").pack() 
    variants = edit_backtracking(matrix, first, second)


def edit_option():
     Label (window, text="Write two sequences: ",font="none 12 bold").pack()
     entry1 = Entry(window, width=150); entry1.pack()   
     entry2 = Entry(window, width=150); entry2.pack() 
     Button(window, text="OK", command =lambda: edit_worker(entry1.get(), entry2.get())).pack()

def clearFrame():
    for widget in window.winfo_children():
        widget.destroy()


def click(value):
   clearFrame()
   if (value==1): fasta_option()
   elif (value == 2): hamming_option()
   elif (value == 3): edit_option()



window = Tk()
window.title("Bioinformatic tools")
window.geometry("800x500+120+120")
Label(window, text="Choose one of these options", font="none 12 bold").pack()

#definitons of Readbuttons for options
r = IntVar()
r1 = Radiobutton(window, text="FASTA file", variable=r, value=1).pack()
r2 = Radiobutton(window, text="Hamming distance", variable=r, value=2).pack()
r3 = Radiobutton(window, text="Sequence alignment using edit distance", variable=r, value=3).pack()
r4 = Radiobutton(window, text="I don't know", variable=r, value=4).pack()
r5 = Radiobutton(window, text="I don't know", variable=r, value=5).pack()
r6 = Radiobutton(window, text="I don't know", variable=r, value=6).pack()
r7 = Radiobutton(window, text="I don't know", variable= r, value=7).pack()

#definition button for ENTER option
enter = Button(window, text="ENTER",command=lambda: click(r.get())).pack()

window.mainloop()