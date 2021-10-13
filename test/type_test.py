import numpy as np

org_file = open('C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\Data\\mycc\\10s\\10s.spe.txt')
print(type(org_file))

table = org_file.read()
length = len(table.split('\n'))
print(length)
entries = table.split('\n')[:-1]

c = []
o = []
for e in entries:
    x = e.split('\t')
    c.append(x[0])
    o.append(x[1])

print(o[0], c[0])
print(o[len(o) - 1], c[len(c) - 1])

contigs_file = open('C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\Data\\mycc\\10s\\10s.fasta')
lines = contigs_file.read().split('>')[1:]
print(len(lines))
print(lines[0])
print(lines[1])
check = np.zeros(len(c))
c = np.array(c, dtype=str)
c.sort()
cont_l = []
for l in lines:
    con = l.split('\n')[0].strip()
    cont_l.append(con)
    arg = c.searchsorted(con)
    check[arg] += 1

cont_arr = np.array(cont_l)
print(cont_arr.shape)
print(np.unique(cont_arr).shape)
print(c.shape)
print(np.unique(c).shape)
print(sum([ce > 1 for ce in check]))
n = sum([ce == 0 for ce in check])

print(length - n)