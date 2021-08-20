import numpy as np

labels = np.load("/home/rom/Dokumente/BA_Data/precomputed/10s/10s_labels.npy")
access = np.load("/home/rom/Dokumente/BA_Data/precomputed/10s/10s_access.npy")
coverage = np.load("/home/rom/Dokumente/BA_Data/precomputed/10s/10s_cov.npy")
counts = np.load("/home/rom/Dokumente/BA_Data/precomputed/10s/10s-count_mat.npy")

print("Access:")
for a in access:
    print(a)

print(len(access))

for l in labels:
    print(l)

print(len(labels))

print("Coverage:")
i = 0
idx = len(coverage)-10
print(idx)
while i < idx:
    print(coverage[i], coverage[i + 1], coverage[i + 2], coverage[i + 3], coverage[i + 4])
    i += 5
print("Coverage_Len:", len(coverage))

print("Counts:")
for c in counts:
    print(len(c))
print("Counts_Len:", len(counts))
print(counts[0])


