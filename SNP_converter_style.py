import re
ol ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}

names = []

for i in range(len(names)):
    s = names[i]
    res = re.split('(\d+)', s)
    if len(res[0]) == 3 and len(res[2]) == 3 and res[2] != 'del' and res[2] != 'dup' :
        print(ol[res[0].upper()]+res[1]+ol[res[2].upper()])    
    elif len(res[0]) == 3 and res[2] == '*':
        print(ol[res[0].upper()]+res[1]+"!")
    elif len(res[0]) == 3 :
        print(ol[res[0].upper()]+res[1]+res[2])
    else:
        print(res[0]+res[1]+res[2])
        

