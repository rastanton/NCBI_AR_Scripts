from Bio import SeqIO
import sys

def OXA_Family_Changer(input_line, family_list):
    List1 = input_line.split('\t')
    Fam = List1[1]
    Allele = List1[0]
    if ('family') in List1[3]:
        Fam = List1[3].split()[0]
        Fam = Fam + '-like'
        if Fam == 'OXA-134-like':
            Fam = 'OXA-235-like'
    elif Allele == 'blaOXA-134':
        Fam = 'OXA-235-like'
    elif (Allele[3:] in family_list):
        Fam = Allele[3:] + '-like'
    return Fam

def OXA_Like_Family_Finder(gene_catalog):
    f = open(gene_catalog, 'r')
    String1 = f.readline()
    List1 = String1[0:-1].split('\t')
    Pos = List1.index('product_name')
    Fams = []
    for line in f:
        if ('blaOXA' in line) and ('family' in line):
            Fam_Line = line.split('\t')[Pos]
            Fam = Fam_Line.split()[0]
            if Fam == 'OXA-134':
                Fam = 'OXA-235'
            if (Fam in Fams) == False:
                Fams.append(Fam)
    return Fams

def Gene_Spreadsheet(input_catalog, output_sheet):
    OXA_Fams = OXA_Like_Family_Finder(input_catalog)
    f = open(input_catalog, 'r')
    Out = open(output_sheet, 'w')
    String1 = f.readline()
    List1 = String1.split('\t')
    Allele_Pos = List1.index('allele')
    Fam_Pos = List1.index('gene_family')
    Subtype_Pos = List1.index('subtype')
    Product_Pos = List1.index('product_name')
    Class_1 = List1.index('class')
    Class_2 = List1.index('subclass')
    for line in f:
        List1 = line.split('\t')
        if List1[Subtype_Pos] == 'AMR':
            Allele = List1[Allele_Pos]
            Fam = List1[Fam_Pos]
            Res_1 = List1[Class_1]
            Res_2 = List1[Class_2]
            if Allele == '':
                Allele = List1[Fam_Pos]
            if List1[Fam_Pos] == 'blaOXA':
                Fam = OXA_Family_Changer(line, OXA_Fams)
                if Fam == 'OXA-134-like':
                    Fam = 'OXA-235-like'
            Out.write(Res_1 + '\t' + Res_2 + '\t' + Fam + '\t' + Allele + '\n')
    f.close()
    Out.close()
                

def Space_Replace(input_line):
    """Replaces a space (' ') with an underscore"""
    Out = input_line
    if (' ' in input_line):
        Out = ''
        for character in input_line:
            if character == ' ':
                Out = Out + '_'
            else:
                Out = Out + character
    return Out
#%%
def AR_DB_Name_Changer(input_fasta, input_catalog, output_fasta):
    """Makes an AR_DB w/ Resistance___sub-resistance__Family__Allele"""
    OXA_Fams = OXA_Like_Family_Finder(input_catalog)
    Lines = []
    f = open(input_catalog, 'r')
    String1 = f.readline()
    for line in f:
        List1 = line.split('\t')
        if List1[6] == 'AMR':
            Lines.append(line)
    f.close()
    Genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Out = open(output_fasta, 'w')
    for gene in Genes:
        Match = gene.description.split()[-1]
        Accession = Match.split(':')[0]
        Pos = Match.split(':')[1]
        Start = Pos.split('-')[0]
        Stop = Pos.split('-')[1]
        Allele = gene.id.split('|')[5]
        Fam = gene.id.split('|')[4]
        for line in Lines:
            Allele = line.split('\t')
            if (Accession in line) and (Start in line) and (Stop in line):
                List1 = line.split('\t')
                Res = List1[7]
                if (' ') in Res:
                    Res = Space_Replace(Res)
                Sub_Res = List1[8]
                if (' ') in Sub_Res:
                    Sub_Res = Space_Replace(Sub_Res)
                if List1[0] == '':
                    Allele = List1[1]
                else:
                    Allele = List1[0]
                Fam = List1[1]
                if ('OXA' in Fam):
                    Fam = OXA_Family_Changer(line, OXA_Fams)
                gene.id = Res + '__' + Sub_Res + '__' + Fam + '__' + Allele
                SeqIO.write(gene, Out, 'fasta')
    Out.close()
    Duplicate_Gene_Renamer(output_fasta, output_fasta)
    CP_Fam_All(output_fasta, output_fasta)

def Total_Counter(input_list, item):
    Count = 0
    for entry in input_list:
        if (entry == item):
            Count += 1
    return Count

def Fasta_Repeat_Renamer(input_fasta, output_fasta):
    """Appends fasta name repeats with numbers"""
    Genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Names = []
    Out_Names = []
    Out = open(output_fasta, 'w')
    for gene in Genes:
        if (gene.id in Names):
            Count = Total_Counter(Names, gene.id)
            Count = str(Count + 1)
            gene.id = gene.id + '_' + Count
            Names.append(gene.id)
        else:
            Names.append(gene.id)
        SeqIO.write(gene, Out, 'fasta')
    Out.close()

def Item_Counter(item, input_list):
    Count = 0
    for items in input_list:
        if item == items:
            Count = Count + 1
    return Count

def Duplicate_Gene_Renamer(input_fasta, output_fasta):
    """Takes in an input fasta and renames any duplicate gene names with a '_count'"""
    Genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Output = open(output_fasta, 'w')
    Name_List = []
    for gene in Genes:
        if (gene.id in Name_List) == True:
            Count = Item_Counter(gene.id, Name_List) + 1
            Name_List.append(gene.id)
            gene.id = gene.id + '_' + str(Count)
            SeqIO.write(gene, Output, 'fasta')
        else:
            Name_List.append(gene.id)
            SeqIO.write(gene, Output, 'fasta')
    Output.close()
        
def Big_5_Fasta(input_AMR_fasta, output_fasta):
    """Makes an AR DB fasta w/ just the big 5 carbapenemases"""
    Genes = list(SeqIO.parse(input_AMR_fasta, 'fasta'))
    Out = open(output_fasta, 'w')
    Bigs = ('VIM', 'IMP', 'KPC', 'NDM', 'OXA-48', 'OXA-235', 'OXA-24', 'OXA-58', 'OXA-23')
    for gene in Genes:
        Fam = gene.id.split('__')[2]
        for entry in Bigs:
            if (entry in Fam):
                SeqIO.write(gene, Out, 'fasta')
    Out.close()

def CP_Fam_Finder(input_fasta):
    Genes = list(SeqIO.parse(input_fasta, 'fasta'))
    CP_Fams = []
    for gene in Genes:
        if ('CARBAPENEM' in gene.id):
            Fam = gene.id.split('__')[-2]
            if (Fam in CP_Fams) == False and Fam != 'blaOXA':
                CP_Fams.append(Fam)
    return CP_Fams

def CP_Fam_All(input_fasta, output_fasta):
##    Fams = CP_Fam_Finder(input_fasta)
    Fams = ['blaVIM', 'blaIMP', 'blaKPC', 'blaNDM', 'OXA-48-like', 'OXA-235-like', 'OXA-24-like', 'OXA-58-like', 'OXA-23-like']
    Genes = list(SeqIO.parse(input_fasta, 'fasta'))
    Out = open(output_fasta, 'w')
    for gene in Genes:
        Fam = gene.id.split('__')[-2]
        if (Fam in Fams) and ('CARBAPENEM' in gene.id) == False:
            Info = gene.id.split('__')
            gene.id = Info[0] + '__CARBAPENEM__' + Info[2] + '__' + Info[3]
        SeqIO.write(gene, Out, 'fasta')
    Out.close()

AR_DB_Name_Changer(sys.argv[1], sys.argv[2], sys.argv[3])
