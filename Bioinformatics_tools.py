from tkinter import *
import math
from tkinter import filedialog
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from Bio.PDB import *
from Bio import pairwise2


#Processing FASTA files
def fasta_option():
    Label(window, text="Give there a FASTA file: ", font="none 12 bold").pack()
    file = filedialog.askopenfilename(initialdir="/", title="Select file")
    sequences = []
    r = IntVar()
    num = 0
    first_dict = {}  # description of sequence
    second_dict = {}  # sequence
    with open(file) as f:
        for line in f:
            if len(line.strip()) != 0:
                splitter = line.split()
                if (splitter[0][:1] == '>'):
                    num += 1
                    sequences.append(splitter[0][1:])
                    string = ' '.join(splitter[1:])
                    first_dict[sequences[num - 1]] = string
                    second_dict[sequences[num - 1]] = None
                else:
                    string = '\n'.join(splitter)
                    if (second_dict[sequences[num - 1]] == None):
                        second_dict[sequences[num - 1]] = string
                    else:
                        second_dict[sequences[num - 1]] = second_dict[sequences[num - 1]] + '\n' + string

    Label(window, text="Choose one of these sequences", font="none 12 bold").pack()
    num = 1
    for sequence in sequences:
        Radiobutton(window, text=sequence, variable=r, value=num).pack()
        num += 1
    Button(window, text="OK", width = 30,command=lambda: fasta_output(num, first_dict, second_dict, sequences, r.get())).pack(anchor = SE)
    Button(window, text="BACK TO MAINMENU", width = 30, command=lambda: main_menu()).pack(anchor = SW)
def fasta_subsequence(name_seq, second_dict, first, second):
    sequence = second_dict.get(name_seq)
    seq = sequence.splitlines()
    sequence = ''.join(seq)
    fir = int(first) - 1
    sec = int(second) - 1
    label = Label(window, text=sequence[fir:sec]).pack()
def fasta_output(num, first_dict, second_dict, sequences, value):
    clearFrame()
    for i in range(0, num):
        if (i == (value - 1)):
            Label(window, text=sequences[i], font="none 12 bold").pack()
            Label(window, text="Description of the sequence", font="Helvetica 10 bold").pack()
            Label(window, text=first_dict.get(sequences[i])).pack()
            Label(window, text="Sequence:", font="Helvetica 10 bold").pack()
            Label(window, text=second_dict.get(sequences[i])).pack()
            Label(window, text="Length of sequence:", font="Helvetica 10 bold").pack()
            Label(window, text=(len(second_dict.get(sequences[i])) + 1)).pack()
            Label(window, text="Subsequence", font="Helvetica 10 bold").pack()
            Label(window, text="Give two numbers between " + '0 and ' + str(len(second_dict.get(sequences[i])) + 1),
                  font="Helvetica 10 bold").pack()
            name = sequences[i]
            entry1 = Entry(window, width=30)
            entry2 = Entry(window, width=30)
            entry1.pack()
            entry2.pack()
            Button(window, text="OK", width = 30, command=lambda: fasta_subsequence(name, second_dict, entry1.get(), entry2.get())).pack(anchor = SE)
            Button(window, text="BACK TO MAINMENU", width = 30, command=lambda: main_menu()).pack(anchor = SW)


#Measuring sequence similarity using Hamming distance
def hamming_worker(first, second): # function which return Hamming distance for two sequences
    res = 0
    if (len(first) == len(second)):
        for i in range(0, len(first)):
            if (first.upper()[i] != second.upper()[i]):
                res += 1
        Label(window, text="Hamming distance between these sequences is: " + str(res), font="none 12 bold").pack()
    else:
        Label(window, text="The two sequences don't have a same length", font="none 12 bold").pack()
def sequences(first, second, file): # function which return two sequences
    first -= 1; second -= 1
    num = 0
    for record in SeqIO.parse(file, "fasta"):
        if (num == first):
            first_seq = record.seq
        elif (num == second):
            second_seq = record.seq
        elif (num == first and num == second):
            first_seq = record.seq
            second_seq = record.seq
        num += 1
    hamming_worker(first_seq, second_seq)
def second_choice(first, file):
    clearFrame()
    Label(window, text="Choose the second sequence", font="none 12 bold").pack()
    r = IntVar()
    num = 1
    for record in SeqIO.parse(file, "fasta"):
        Radiobutton(window, text=record.id, variable=r, value=num).pack()
        num += 1
    Button(window, text="OK", width=30, command=lambda: sequences(first,r.get(), file)).pack(anchor=SE)
    Button(window, text="BACK", width=30, command=lambda: hamming_option()).pack(anchor=SW)
def hamming_choice(value):
    if (value == 1): # input is a FASTA format file
        clearFrame()
        Label(window, text="Give there a FASTA file: ", font="none 12 bold").pack()
        file = filedialog.askopenfilename(initialdir="/", title="Select file")
        clearFrame()
        Label(window, text="Choose the first sequence", font="none 12 bold").pack()
        r = IntVar()
        num = 1
        for record in SeqIO.parse(file, "fasta"):
            Radiobutton(window, text=record.id, variable=r, value=num).pack()
            num += 1
        Button(window, text="OK", width=30, command=lambda: second_choice(r.get(), file)).pack(anchor=SE)
        Button(window, text="BACK", width=30, command=lambda: hamming_option()).pack(anchor=SW)
    else: # input is manually entered sequence
        clearFrame()
        Label(window, text="Write the first sequence:", font="none 12 bold").pack()
        entry1 = Entry(window, width=150);
        entry1.pack()
        Label(window, text="Write the second sequence:", font="none 12 bold").pack()
        entry2 = Entry(window, width=150);
        entry2.pack()
        Button(window, text="OK", width=30, command=lambda: hamming_worker(entry1.get(), entry2.get())).pack(anchor = SE)
        Button(window, text="BACK TO MAINMENU", width=30, command=lambda: main_menu()).pack(anchor = SW)
def hamming_option(): # function which receive an input
    clearFrame()
    Label(window, text="Choose if you want to load a FASTA file or only write two sequences:",
          font="none 12 bold").pack()
    r = IntVar()
    r1 = Radiobutton(window, text="FASTA file", variable=r, value=1).pack()
    r2 = Radiobutton(window, text="Write two sequences", variable=r, value=2).pack()
    Button(window, text="OK", width=30,command=lambda: hamming_choice(r.get())).pack(anchor = SE)
    Button(window, text="BACK TO MAINMENU", width = 30, command=lambda: main_menu()).pack(anchor = SW)

#Sequence alignment using edit distance
def edit_backtracking(matrix, first, second): # function which return and variants of alignment
    text = make_textbox()
    Label(window, text="Sequence alignment: ", font="none 12 bold").pack()
    for a in pairwise2.align.globalxx(first, second):
        text.insert(END, format_alignment(*a))
def make_matrix(first, second):
    matrix = []
    for i in range(first):
        matrix.append([0] * second)
    return matrix
def edit_worker(first, second):
    matrix = make_matrix(len(first) + 1, len(second) + 1)
    num = 0
    # first row
    for i in range(len(second) + 1):
        matrix[0][i] = num
        num += 1
    num = 0
    # first column
    for j in range(len(first) + 1):
        matrix[j][0] = num
        num += 1
    for m in range(1, len(first) + 1):
        for n in range(1, len(second) + 1):
            if (first.upper()[m - 1] == second.upper()[n - 1]):
                matrix[m][n] = matrix[m - 1][n - 1]
            else:
                matrix[m][n] = min(matrix[m - 1][n - 1] + 1, matrix[m][n - 1] + 1, matrix[m - 1][n] + 1)
        edit_distance = matrix[m][n]
    Label(window, text="Edit distance: " + str(edit_distance), font="none 12 bold").pack()
    variants = edit_backtracking(matrix, first, second)
    Button(window, text="BACK TO MAINMENU", width = 30, command=lambda: main_menu()).pack(anchor = SW)
def edit_option(): # funtion for receive an input
    Label(window, text="Write two sequences: ", font="none 12 bold").pack()
    entry1 = Entry(window, width=150);
    entry1.pack()
    entry2 = Entry(window, width=150);
    entry2.pack()
    Button(window, text="OK", width = 30, command=lambda: edit_worker(entry1.get(), entry2.get())).pack(anchor = SE)
    Button(window, text="BACK TO MAINMENU", width = 30, command=lambda: main_menu()).pack(anchor = SW)

#Processing PDB files
def pdb_option(): # print all PDB options
    clearFrame()
    Label(window, text="Give there a PDB file: ", font="none 12 bold").pack()
    file = filedialog.askopenfilename(initialdir="/", title="Select file")
    r = IntVar()
    r1 = Radiobutton(window, text="Object reprezenting a model", variable=r, value=1).pack()
    r2 = Radiobutton(window, text="Object reprezenting a structure(chain) within a model", variable=r, value=2).pack()
    r3 = Radiobutton(window, text="Object reprezenting residuum within a chain", variable=r, value=3).pack()
    r4 = Radiobutton(window, text="Object reprezenting an atom within a residue", variable=r, value=4).pack()
    r5 = Radiobutton(window, text="Information about the stored structure", variable=r, value=5).pack()
    r6 = Radiobutton(window, text="The width of the structure", variable=r, value=6).pack()
    r7 = Radiobutton(window, text="List of atoms being in given distance from given ligand", variable=r, value=7).pack()
    enter = Button(window, text="OK", width = 30, command=lambda: click_pdb(r.get(), file)).pack(anchor = SE)
    Button(window, text="BACK TO MAINMENU", width=30, command=lambda: main_menu()).pack(anchor = SW)
def pdb_parser(file):
    atoms = {} # where keys = atom serial number and atom's name and values = x,y,z
    hetatm = {} # where keys = hetatm serial number and hetatm name and values = x,y,z

    with open (file) as f:
        for line in f:
            splitter = line.split()
            if (splitter[0] == 'ATOM'):
                atoms[splitter[1]] = splitter[6:9]
                val = atoms.get(splitter[1])
                val.append(splitter[2])
                atoms[splitter[1]] = val
            elif(splitter[0] == 'HETATM'):
                hetatm[splitter[1]] = splitter[6:9]
                val = hetatm.get(splitter[1])
                val.append((splitter[2]))
                hetatm[splitter[1]] = val
    return (atoms, hetatm)
def give_distance(atom_coordinates, hetatm_coordinates):
    a = (float(atom_coordinates[0]) - float(hetatm_coordinates[0]))**2
    b = (float(atom_coordinates[1]) - float(hetatm_coordinates[1]))**2
    c = (float(atom_coordinates[2]) - float(hetatm_coordinates[2]))**2
    return (math.sqrt(a+b+c))
def give_atoms(atoms, hetatm, val, num):
    clearFrame()
    atom_right = [] # atoms which are in given distance from given ligand
    hetatm_coordinates = hetatm.get(val)
    hetatm_coordinates = hetatm_coordinates[0:3]
    #print(hetatm_coordinates)
    for keys in atoms:
        atom_coordinates = atoms.get(keys)[0:3]
        distance = give_distance(atom_coordinates, hetatm_coordinates)
        if (distance <= float(num)):
            atom_right.append(keys)
    Button(window, text = "BACK", command=lambda: pdb_option()).pack()
    Label(window, text = "Atoms being in given distance from given ligand: ", font="none 12 bold").pack()
    text = Text(window, width=100)
    text.pack(side=LEFT, fill=Y)
    scroll = Scrollbar(window)
    scroll.pack(side=RIGHT, fill=Y)
    scroll.config(command=text.yview())
    text.config(yscrollcommand=scroll.set)
    for a in atom_right:
        value = atoms.get(a)
        value = value[3]
        text.insert(END, str(a) + " " + str(atoms.get(a)[3]) + "\n")
def get_model(file): # function which return list of models from PDB
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = parser.get_structure("test", file)
    model = structure.get_models()
    models = list(model)
    return models
def get_max_distance(atoms):
    max_distance = 0
    coordinates = []
    center = [] #coordinates of center
    for i in atoms:
        atom_coordinates = atoms.get(i)[0:3]
        for j in atoms:
            atom_coordinates2 = atoms.get(j)[0:3]
            distance = give_distance(atom_coordinates, atom_coordinates2)
            if (distance > max_distance):
                max_distance = distance
                x = abs(float(atom_coordinates[0]) - float(atom_coordinates2[0])) / 2
                y = abs(float(atom_coordinates[1])- float(atom_coordinates2[1])) / 2
                z = abs(float(atom_coordinates[2]) - float(atom_coordinates2[2])) / 2
                center.append(x); center.append(y); center.append(z)
                coordinates.clear()
                coordinates.append(atom_coordinates); coordinates.append(atom_coordinates2)
    return (max_distance, center, coordinates)

def click_pdb(value, file): # PDB parser
    clearFrame()
    models = get_model(file)
    models_len = len(models)
    chains_len = 0; residues_len = 0; atoms_len = 0
    pdb_parser(file)
    for m in models:
        chains = list(m.get_chains())
        chains_len += len(chains)
        for ch in chains:
            residue = list(ch.get_residues())
            residues_len += len(residue)
            for r in residue:
                atoms = list(r.get_atoms())
                atoms_len += len(atoms)
    if (value == 1): # Write the models of PDB file
        clearFrame()
        Label(window, text="Object representing a model", font="none 12 bold").pack()
        for m in models:
            Label(window, text=m).pack()
        Button(window, text="BACK", width = 30, command=lambda: pdb_option()).pack()
    elif (value == 2): # Write chains of the PDB file
        clearFrame()
        Label(window, text="Object reprezenting a structure(chain) within a model", font="none 12 bold").pack()
        for m in models:
            chains = list(m.get_chains())
            chains_len = len(chains)
            for ch in chains:
                Label(window, text=str(m)+" "+str(ch)).pack()
        Button(window, text="BACK", width = 30, command=lambda: pdb_option()).pack()
    elif(value == 3): # Write the residuums in PDB file
        clearFrame()
        Label(window, text="Object reprezenting residuum within a chain", font="none 12 bold").pack()
        Button(window, text="BACK", width=30, command=lambda: pdb_option()).pack()
        text = make_textbox()
        for m in models:
            chains = list(m.get_chains())
            for ch in chains:
                residue = list(ch.get_residues())
                residuums_len = len(residue)
                for r in residue:
                    text.insert(END, str(ch) + " " + str(r) + "\n")
    elif(value == 4): # Write the atoms in PDB file
        clearFrame()
        Label(window, text="Object reprezenting an atom within a residue", font="none 12 bold").pack()
        Button(window, text="BACK", width=30, command=lambda: pdb_option()).pack()
        text = make_textbox()
        for m in models:
            chains = list(m.get_chains())
            for ch in chains:
                residue = list(ch.get_residues())
                for r in residue:
                    atoms = list(r.get_atoms())
                    atoms_len = len(atoms)
                    for a in atoms:
                        text.insert(END, str(r) + " " + str(a) + "\n")
    elif(value == 5): # Write summary information about PDB file
        clearFrame()
        Label(window, text="Information about the stored structure", font="none 12 bold").pack()
        Label(window, text="Number of models: " + str(models_len)).pack()
        Label(window, text="Number of structures: " + str(chains_len)).pack()
        Label(window, text="Number of residues: " + str(residues_len)).pack()
        Label(window, text="Number of atoms: " + str(atoms_len)).pack()
        Button(window, text="BACK", command=lambda: pdb_option()).pack()
    elif (value == 6): # Write the width of the structure in PDB
        clearFrame()
        Label(window, text="The width of the structure", font="none 12 bold").pack()
        atoms, hetatm = pdb_parser(file)
        distance, center, coordinates = get_max_distance(atoms)
        Label(window, text=distance).pack()


        Button(window, text="BACK", width=30, command=lambda: pdb_option()).pack()
    elif (value == 7): # Write the HETATM atoms
        atoms, hetatm = pdb_parser(file)
        Label(window, text="Write number reprezenting distance from choosen ligand", font="none 12 bold").pack()
        entry1 = Entry(window, width=30)
        entry1.pack()
        Label(window, text="Write one of these number of ligands: ", font="none 12 bold").pack()
        entry2 = Entry(window, width=30)
        entry2.pack()
        Button(window, text="OK", width=30, command=lambda: give_atoms(atoms, hetatm, entry2.get(), entry1.get())).pack(
            anchor=SE)
        Button(window, text="BACK", width=30, command=lambda: pdb_option()).pack(anchor=SW)
        text = make_textbox()
        for key in hetatm:
            name = str(key) + " " + str(hetatm.get(key)[3])
            text.insert(END, str(name) + '\n')


#Processing multiple sequence alignment
def get_scoring_matrix(scoring_file): # function which process the scoring matrix
    scoring_matrix = []
    with open(scoring_file) as f:
        for line in f:
            num = -1
            splitter = line.split(';')
            scoring_matrix.append(splitter[:len(splitter) - 1])
    return (scoring_matrix)
def get_scores (scoring_mat, adict): # function which return score for substitution in MSA
    analphabets = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    column_scores = []
    res = 0
    gap = 0.5
    for value in adict.values():
        length = len(value[0])  # length of block
        length2 = len(value)  # number of sequences
        for i in range(0, length2):
            col_score = 0
            for j in range(0, length):
                for k in range (j+1, length):
                    if (value[i][j] != '-' and value[i][k] != '-'):
                        index1 = analphabets.index(value[i][j])
                        index2 = analphabets.index(value[i][k])
                        col_score += int(scoring_mat[index1][index2])
            column_scores.append(col_score)
    return(column_scores)
def get_dictionary(file):  #work with the CLUSTAL file format
    block = []  # one Clustal block
    adict = {} # dictionary with Clustal MSA blocks
    num = 0
    with open(file) as f:
        for line in f:
            if (len(line.strip()) != 0 and line[0] != "C"):
                splitter = line.split()
                if (splitter[0][0] != '*' and splitter[0][0] != '.' and splitter[0][0] != ':'):
                    seq = splitter[1]
                    if (len(block) == 0):
                        for i in range(0, len(seq)):
                            block.append([seq[i]])
                    else:
                        for i in range(0, len(seq)):
                            block[i].append(seq[i])
            if line in ['\n', '\r\n']:
                if not block:
                    continue
                else:
                    adict[num] = block
                    block = []
                    num += 1
        adict[num] = block
        length = len(block)
    return (adict, length)
def get_sequence_worker(value,sequence, adict, file): # function which print the right sequence
    clearFrame()
    button = Button(window, text="BACK", command=lambda: get_sequence(value, adict, file)).pack(anchor=SW)
    text = make_textbox()
    if (sequence.isnumeric()):
        if (int(sequence) <= len(adict.keys())):
            key = list(adict.keys())[int(sequence) - 1]
            Label(window, text="Sequence: ", font="none 12 bold").pack()
            text.insert(END, adict.get(key))
        else:
            Label(window, text="Wrong number of sequence.", font="none 12 bold").pack()
    else:
        if (adict[sequence] == None):
            Label(window, text="Wrong name of sequence.", font="none 12 bold").pack()
        else:
            Label(window, text="Sequence: " + adict.get(sequence), font="none 12 bold").pack()
def get_sequence(value, adict, file): # function which recieve users input
    clearFrame()
    entry = Entry(window, width=150);
    Label(window, text="Write sequence position or ID", font="none 12 bold").pack()
    entry.pack()
    Button(window, text="OK", width = 30, command=lambda: get_sequence_worker(value, entry.get(), adict, file)).pack(anchor = SE)
    button = Button(window, text="BACK", width=30, command=lambda: multiple_seq_option()).pack(anchor = SW)
def get_column_worker(num, adict, file): # function which print the right column from MSA
    clearFrame()
    key = list(adict.keys())[0] # one column from MSA
    length = len(adict.get(key)) # length of one column in MSA
    Button(window, text="BACK", width=30, command=lambda: get_column(adict, file)).pack()
    text = make_textbox()
    if (int(num) < 0 and int(num) >= length):
        Label(window, text="Wrong number.", font="none 12 bold").pack()
    else:
        Label(window, text="The " + num + "column is: ", font="none 12 bold").pack()
        adict2, length = get_dictionary(file) # adict2 - dictionary with MSA blocks, length - length of one column in MSA
        column = 0
        num = int(num); num -= 1
        while(num > length):
            num -= length
            column += 1
        text.insert(END, str(adict2.get(column)[num]) + '\n')
def get_column(adict, file): # function which only receive an input
    clearFrame()
    entry = Entry(window, width=150)
    Label(window, text="Write number of column of MSA", font="none 12 bold").pack()
    entry.pack()
    Button(window, text="OK", width = 30, command=lambda: get_column_worker(entry.get(), adict, file)).pack(anchor = SE)
    Button(window, text="BACK", width = 30, command=lambda: multiple_seq_option()).pack(anchor = SW)
def get_sum_of_pairs(file): #fuction for evaluating sum of pairs of MSA
    clearFrame()
    res = 0
    lab = Label(window, text = "Load the scoring matrix", font="none 12 bold").pack()
    scoring_file = filedialog.askopenfilename(initialdir="/", title="Select file")
    scoring_matrix = get_scoring_matrix(scoring_file)
    adict, length = get_dictionary(file) # adict - dictionary with MSA blocks, length - length of one column in MSA
    column_scores = get_scores(scoring_matrix, adict)
    for score in column_scores:
        res += score
    Label(window, text="Score of MSA is: " + str(res), font="none 12 bold").pack()
    button = Button(window, text="BACK", width=30, command=lambda: multiple_seq_option()).pack()
def click3(value, adict, file): #function with 3 options
    clearFrame()
    if (value == 1):
        get_sequence(value, adict, file)
    elif (value == 2):
        get_column(adict, file)
    elif (value == 3):
        get_sum_of_pairs(file)
def multiple_seq_option():
    clearFrame()
    sequences = [] # names of sequences
    adict = {} # dictionary where keys are names of sequences and values are sequences
    Label(window, text="Give a CLUSTAL file: ", font="none 12 bold").pack()
    file = filedialog.askopenfilename(initialdir="/", title="Select file")
    with open(file) as f:
        for line in f:
            if (line[0] == "U"):
                splitter = line.split()
                sequences.append(splitter[0])
                if (splitter[0] in adict):
                    adict[splitter[0]] = adict[splitter[0]] + '\n' + splitter[1]
                else:
                    adict[splitter[0]] = splitter[1]
    Label(window, text="Choose one of these options", font="none 12 bold").pack()
    r = IntVar()
    r1 = Radiobutton(window, text="Retrieve sequence by its position or ID", variable=r, value=1).pack()
    r2 = Radiobutton(window, text="Retrieve given column from the MSA", variable=r, value=2).pack()
    r3 = Radiobutton(window, text="Retrieve sum of pairs score of a column and whole MSA with respect to given scoring matrix", variable=r, value=3).pack()
    enter = Button(window, text="OK", width = 30, command=lambda: click3(r.get(), adict, file)).pack(anchor=SE)
    Button(window, text="BACK TO MAINMENU", width = 30, command=lambda: main_menu()).pack(anchor = SW)

#Conservation determination from MSA
def get_consensus(adict):
    consensus = []
    counter = {}
    for key in adict:
        seq = adict.get(key)
        for s in seq:
            length = len(s)
            counter.clear()
            num = 0
            for n in range (0, length):
                if (s[n] in counter.keys()):
                    val = counter.get(s[n])
                    counter[s[n]] = val + 1
                else:
                    counter[s[n]] = 1
            for c in counter:
                if (counter.get(c) > (length/2)):
                    consensus.append(c)
                    num = 1
            if (num == 0):
                consensus.append('-')

    return(consensus)
def conservation_option():
    Label(window, text="Give a CLUSTAL file: ", font="none 12 bold").pack()
    file = filedialog.askopenfilename(initialdir="/", title="Select file")
    num = 0
    adict, length = get_dictionary(file)
    Button(window, text="BACK TO MAIN MENU", width=30, command=lambda: main_menu()).pack(anchor=SW)
    Label(window, text="Give a Scoring matrix file: ", font="none 12 bold").pack()
    file2 = filedialog.askopenfilename(initialdir="/", title="Select file")
    scoring_matrix = get_scoring_matrix(file2)
    column_scores = get_scores(scoring_matrix, adict)
    max_score = 0
    indexes_top_columns = []
    for i in column_scores:
        if (i >= max_score):
            max_score = i
    for j in column_scores:
        if (j == max_score):
            indexes_top_columns.append(column_scores.index(j))
    consensus = get_consensus(adict)
    consensus_string = ''.join(consensus)
    Label(window, text="Consensus sequence and the best scoring columns are: ", font="none 12 bold").pack()
    text = make_textbox()
    text.insert(END, "Consensus sequence: " + '\n')
    text.insert(END, str(consensus_string) + '\n')
    text.insert(END, "\n" +"Indexes of top scoring columns: ")
    for i in indexes_top_columns:
        text.insert(END, str(i+1) + " ")
    text.insert(END, "\n")
    num = 0
    for j in indexes_top_columns:
        j-= 1
        while(j > length):
            j -= length
            num += 1
        text.insert(END, str(adict.get(num)[j]) + '\n')


#Compute strucure-related properties
def constants_maker(a, b, c):
    return(a**2*(b**2 + c**2 - a**2))
def get_center(x1,x2,x3,alfa, beta,gamma):
    return((alfa*x1 + beta*x2 + gamma*x3)/(alfa + beta + gamma))
def give_diameter(file):
    clearFrame()
    atoms, hetatms = pdb_parser(file)
    diameter = 0
    distance, center, coordinates = get_max_distance(atoms)
    radius = distance / 2
    for i in atoms:
        atom_coordinates = atoms.get(i)[0:3]
        distance2 = give_distance(center, atom_coordinates)
        if (distance2 > distance):
            # I have 3 points and I want have a circumcenter in 3D
            a = give_distance(coordinates[1], atom_coordinates)
            b = give_distance(atom_coordinates, coordinates[0])
            c = give_distance(coordinates[0], coordinates[1])
            alfa = constants_maker(a, b, c)
            beta = constants_maker(b, c, a)
            gamma = constants_maker(c, b, a)
            center_x = (coordinates[0][0], coordinates[1][0], atom_coordinates[0], alfa, beta, gamma)
            center_y = (coordinates[0][1], coordinates[1][1], atom_coordinates[1], alfa, beta, gamma)
            center_z = (coordinates[0][2], coordinates[1][2], atom_coordinates[2], alfa, beta, gamma)
            radius *= 1.0001
            if (diameter < (radius * 2)):
                diameter = radius * 2
    Label(window, text="Diameter of this protein is: " + str(diameter)).pack()
    Button(window, text="BACK", width=30, command=lambda: structure_option()).pack(anchor=SW)
def get_ratio(file):
    neighbors = {} # where key = atom and value = number of neighbours
    structure = PDBParser().get_structure('test', file)
    atoms = Selection.unfold_entities(structure, 'A')
    buried_atoms = []
    surface_atoms = []
    ns = NeighborSearch(atoms)
    backbones_atoms = ["N", "CA", "C", "O"]
    for i in range (0, len(atoms)):
        center = atoms[i].get_coord()
        num_neighbors = len(ns.search(center, 2.75)) #distance 2.5 A, because it is van der wals radius for water
        if (num_neighbors < 5):
            if (atoms[i].name in backbones_atoms):
                buried_atoms.append(atoms[i])
            else:
                surface_atoms.append(atoms[i])
        else:
            buried_atoms.append(atoms[i])
    num_buried = len(buried_atoms)
    num_surface = len(surface_atoms)
    return(buried_atoms, surface_atoms)

def give_statistic (adict):
    atoms = {}
    for key in adict:
        if key.name in atoms:
            atoms[key.name] = atoms.get(key.name) + 1
        else:
            atoms[key.name] = 1
    return atoms

def sum_of_amk(adict):
    num = 0
    for key in adict:
        num += adict.get(key)
    return (num)

def get_polar_amk(adict):
    polar_AMK = ["S", "T", "C", "D", "G", "N", "E", "K", "R", "H"]
    polar_amk = {}  # where key = name of atom and value = sum of these atoms
    num_amk = 0
    for key in adict:
        if key in polar_AMK:
            num_amk += adict.get(key)
            polar_amk[key] = adict.get(key)
    return(polar_amk, num_amk)


def structure_click(val, file):
    clearFrame()
    buried_atoms, surface_atoms = get_ratio(file)
    buried_stat = give_statistic(buried_atoms)  # statistics of buried atoms
    surface_stat = give_statistic(surface_atoms)  # statistics of surface atoms
    if (val == 1):
        give_diameter(file)
    elif(val == 2):
        Label(window, text="The ratio of surface and buried amino acids: ").pack()
        Label(window, text="buried: " + str(len(buried_atoms))).pack()
        Label(window, text="surface: " + str(len(surface_atoms))).pack()
        Label(window, text="ratio: " + str(len(buried_atoms) / len(surface_atoms))).pack()
        Button(window, text="BACK", width=30, command=lambda: structure_option()).pack(anchor=SW)
    elif(val == 3):
        Label(window, text="Statistic information for histogram: ").pack()
        Button(window, text="BACK", width=30, command=lambda: structure_option()).pack(anchor=SW)
        text = make_textbox()
        text.insert(END, "Buried amino acids: " + "\n")
        for atoms in buried_stat:
            text.insert(END, str(atoms) + " " + str(buried_stat.get(atoms)) + ";")
        text.insert(END, "\n" + "Surface amino acids: " + "\n")
        for atoms in surface_stat:
            text.insert(END, str(atoms) + " " + str(surface_stat.get(atoms)) + ";")
    elif(val == 4):
        polar_buried, num_buried = get_polar_amk(buried_stat)
        polar_surface, num_surface = get_polar_amk(surface_stat)
        Label(window, text="Portion of polar aminoacids: ").pack()
        Button(window, text="BACK", width=30, command=lambda: structure_option()).pack(anchor=SW)
        text = make_textbox()
        text.insert(END, "Polar surface amino acids: " + str(num_surface) + "\n")
        for i in polar_surface:
            text.insert(END, str(i) + ":" + str(polar_surface.get(i)) + " ;")
        text.insert(END, "\n" + "Buried surface amino acids: " + str(num_buried) + "\n")
        for i in polar_buried:
            text.insert(END, str(i) + ":" + str(polar_buried.get(i)) + " ;")
    elif(val == 5):
        Label(window, text="Give a PDB file: ").pack()
        file2 = filedialog.askopenfilename(initialdir="/", title="Select file")
        clearFrame()
        buried_atoms2, surface_atoms2 = get_ratio(file2)
        buried_stat2 = give_statistic(buried_atoms2)
        surface_stat2 = give_statistic(surface_atoms2)
        num_buried = sum_of_amk(buried_stat); num_buried2 = sum_of_amk(buried_stat2)
        num_surface = sum_of_amk(surface_stat); num_surface2 = sum_of_amk(surface_stat2)
        Button(window, text="BACK", width=30, command=lambda: structure_option()).pack(anchor=SW)
        Label(window, text="Buried and surface amino acids: " + "\n").pack()
        Label(window, text="First PDB file: " + "\n").pack()
        Label(window, text="Buried: " + str(num_buried) + "\n").pack()
        Label(window, text="Surface: " + str(num_surface) + "\n").pack()
        Label(window, text="Second PDB file: " + "\n").pack()
        Label(window, text="Buried: " + str(num_buried2) + "\n").pack()
        Label(window, text="Surface: " + str(num_surface2) + "\n" + "\n").pack()
        polar_buried, num_pburied = get_polar_amk(buried_stat)
        polar_surface, num_psurface = get_polar_amk(surface_stat)
        polar_buried2, num_pburied2 = get_polar_amk(buried_stat2)
        polar_surface2, num_psurface2 = get_polar_amk(surface_stat2)
        Label(window, text="Polar buried and surface amino acids: " + "\n").pack()
        text = make_textbox()
        text.insert(END, "The first PDB file: " + "\n")
        text.insert(END, "Polar surface amino acids: " + str(num_psurface) + "\n")
        for i in polar_surface:
            text.insert(END, str(i) + ":" + str(polar_surface.get(i)) + " ;")
        text.insert(END, "\n" + "Buried surface amino acids: " + str(num_pburied) + "\n")
        for i in polar_buried:
            text.insert(END, str(i) + ":" + str(polar_buried.get(i)) + " ;")
        text.insert(END, "\n" + "\n" + "The second PDB file: " + "\n")
        text.insert(END, "Polar surface amino acids: " + str(num_psurface2) + "\n")
        for i in polar_surface2:
            text.insert(END, str(i) + ":" + str(polar_surface2.get(i)) + " ;")
        text.insert(END, "\n" + "Buried surface amino acids: " + str(num_pburied2) + "\n")
        for i in polar_buried2:
            text.insert(END, str(i) + ":" + str(polar_buried2.get(i)) + " ;")

def structure_option():
    clearFrame()
    Label(window, text="Give a PDB file: ").pack()
    file = filedialog.askopenfilename(initialdir="/", title="Select file")
    r = IntVar()
    r1 = Radiobutton(window, text="Compute the diameter of the protein", variable=r, value=1).pack()
    r2 = Radiobutton(window, text="Compute the ratio of surface and buried amino acids", variable=r, value=2).pack()
    r3 = Radiobutton(window, text="Output histogram of amino acids composition of buried and exposed amino acids", variable=r, value=3).pack()
    r4 = Radiobutton(window, text="Quantify portion of polar amino acids in the core and on the surface od the protein", variable=r, value=4).pack()
    r5 = Radiobutton(window, text="Use the second structure to quantify the portion surface and buried amino acids", variable=r, value=5).pack()
    enter = Button(window, text="OK", width=30, command=lambda: structure_click(r.get(), file)).pack(anchor=SE)
    button = Button(window, text="BACK TO MAINMENU", width=30, command=lambda: main_menu()).pack(anchor=SW)

#main function for GUI and the main window
def clearFrame():
    for widget in window.winfo_children():
        widget.destroy()
def make_textbox():
    text = Text(window, width=100)
    text.pack(side=LEFT, fill=Y)
    scroll = Scrollbar(window)
    scroll.pack(side=RIGHT, fill=Y)
    scroll.config(command=text.yview())
    text.config(yscrollcommand=scroll.set)
    return(text)

def click(value):
    clearFrame()
    if (value == 1):
        fasta_option()
    elif (value == 2):
        hamming_option()
    elif (value == 3):
        edit_option()
    elif(value == 4):
        pdb_option()
    elif(value == 5):
        multiple_seq_option()
    elif(value == 6):
        conservation_option()
    elif(value == 7):
        structure_option()
def main_menu():
    clearFrame()
    Label(window, text="Choose one of these options", font="none 12 bold").pack()
    r = IntVar()
    r1 = Radiobutton(window, text="FASTA parser", variable=r, value=1).pack()
    r2 = Radiobutton(window, text="Hamming distance", variable=r, value=2).pack()
    r3 = Radiobutton(window, text="Sequence alignment using edit distance", variable=r, value=3).pack()
    r4 = Radiobutton(window, text="PDB parser", variable=r, value=4).pack()
    r5 = Radiobutton(window, text="Processing multiple sequnce alignment", variable=r, value=5).pack()
    r6 = Radiobutton(window, text="Conservation determination from MSA", variable=r, value=6).pack()
    r7 = Radiobutton(window, text="Compute structure-related properties", variable=r, value=7).pack()
    enter = Button(window, text="ENTER", command=lambda: click(r.get())).pack()



window = Tk()
window.title("Bioinformatics tools")
window.geometry("800x500+120+120")
main_menu()
window.mainloop()