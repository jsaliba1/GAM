strain = '' #+ or -
rg = '' #DNA Seq
sm = 0#start pos
em = 0#edn pos

pos_list = [] #pos list 
ref_list = [] #DNA ref
alt_list = [] #DNA alt


rg = list(rg.upper())
refs = list(zip(*[iter(rg)]*3))

amino = {
    'TCA': 'S',    # Serina
    'TCC': 'S',    # Serina
    'TCG': 'S',    # Serina
    'TCT': 'S',    # Serina
    'TTC': 'F',    # Fenilalanina
    'TTT': 'F',    # Fenilalanina
    'TTA': 'L',    # Leucina
    'TTG': 'L',    # Leucina
    'TAC': 'Y',    # Tirosina
    'TAT': 'Y',    # Tirosina
    'TAA': '!',    # Stop
    'TAG': '!',    # Stop
    'TGC': 'C',    # Cisteina
    'TGT': 'C',    # Cisteina
    'TGA': '!',    # Stop
    'TGG': 'W',    # Triptofano
    'CTA': 'L',    # Leucina
    'CTC': 'L',    # Leucina
    'CTG': 'L',    # Leucina
    'CTT': 'L',    # Leucina
    'CCA': 'P',    # Prolina
    'CCC': 'P',    # Prolina
    'CCG': 'P',    # Prolina
    'CCT': 'P',    # Prolina
    'CAC': 'H',    # Histidina
    'CAT': 'H',    # Histidina
    'CAA': 'Q',    # Glutamina
    'CAG': 'Q',    # Glutamina
    'CGA': 'R',    # Arginina
    'CGC': 'R',    # Arginina
    'CGG': 'R',    # Arginina
    'CGT': 'R',    # Arginina
    'ATA': 'I',    # Isoleucina
    'ATC': 'I',    # Isoleucina
    'ATT': 'I',    # Isoleucina
    'ATG': 'M',    # Methionina
    'ACA': 'T',    # Treonina
    'ACC': 'T',    # Treonina
    'ACG': 'T',    # Treonina
    'ACT': 'T',    # Treonina
    'AAC': 'N',    # Asparagina
    'AAT': 'N',    # Asparagina
    'AAA': 'K',    # Lisina
    'AAG': 'K',    # Lisina
    'AGC': 'S',    # Serina
    'AGT': 'S',    # Serina
    'AGA': 'R',    # Arginina
    'AGG': 'R',    # Arginina
    'GTA': 'V',    # Valina
    'GTC': 'V',    # Valina
    'GTG': 'V',    # Valina
    'GTT': 'V',    # Valina
    'GCA': 'A',    # Alanina
    'GCC': 'A',    # Alanina
    'GCG': 'A',    # Alanina
    'GCT': 'A',    # Alanina
    'GAC': 'D',    # Acido Aspartico
    'GAT': 'D',    # Acido Aspartico
    'GAA': 'E',    # Acido Glutamico
    'GAG': 'E',    # Acido Glutamico
    'GGA': 'G',    # Glicina
    'GGC': 'G',    # Glicina
    'GGG': 'G',    # Glicina
    'GGT': 'G'     # Glicina
}


if strain == '+':
    for i in range(len(pos_list)):
        if rg[pos_list[i]-sm] == ref_list[i]:
            nrg = rg.copy()
            nrg[pos_list[i]-sm] = alt_list[i]
            
            alts = list(zip(*[iter(nrg)]*3))
            
            for j in range(len(refs)):
                if refs[j] != alts[j]:
                    print(amino[''.join(refs[j])]+str(j+1)+amino[''.join(alts[j])])

if strain == '-':
    for i in range(len(pos_list)):
        if rg[em - pos_list[i]] == ref_list[i]:
            nrg = rg.copy()
            nrg[em- pos_list[i]] = alt_list[i]
            
            alts = list(zip(*[iter(nrg)]*3))
            
            for j in range(len(refs)):
                if refs[j] != alts[j]:
                    print(amino[''.join(refs[j])]+str(j+1)+amino[''.join(alts[j])])