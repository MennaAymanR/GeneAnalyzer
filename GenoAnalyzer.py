
def transcription(user_input):
    for i in user_input:
        if i == "A":
            print("U")
        if i=="T":
            print("A")
        if i=="C":
            print("G")
        if i=="G":
            print("C")
        if i == ("N"):
            print("AnyNucloetide")
def translate(user_input):
    dictCodon = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M','ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T','AAC': 'N',
        'AAT': 'N', 'AAA': 'K', 'AAG': 'K','AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L','CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q','CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V','GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E','GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S','TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_','TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    protein = " "
    if len(user_input) % 3 == 0:
        for i in range(0, len(user_input), 3):
            codon = user_input[i:i + 3]
            protein += dictCodon[codon]
        return protein
def filter(user_input):

    A=user_input.count("A")
    C=user_input.count("C")
    G=user_input.count("G")
    T=user_input.count("T")
    N=user_input.count("N")
    print("the ratio of N bases is:",N/(N+A+C+G+T) * 100)
    FilteredSeq = user_input.replace("N" , "")
    print(FilteredSeq)
def start_code(user_input):
    for i in range(len(user_input)):
        if user_input[i]=="A":
            if user_input[i+1]=="T":
                if user_input[i + 2] == "G":
                    print(i)
def eng_name():
    print("developed by menna ayman ")

user_input = input("Do you input a path or seq? ")
seq = ""
if user_input == "path":
    x = input("Please enter file's path")
    y = open(x,"r")
    seq = y.read()
    print(seq)
elif user_input == "seq":
    seq = input("please enter the sequence:")
    seq = seq.upper()
    print(seq)
X = input("what's your process?transcripe or translate or filter or start_code or eng_name: ")
if X == "transcripe":
    print(X)
    transcription(seq)
elif X == "translate":
    translate(seq)
elif X == "filter":
    filter(seq)
elif X == "start_code":
    start_code(seq)
elif x == "eng_name":
    eng_name()


