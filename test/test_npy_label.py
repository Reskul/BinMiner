import numpy as np

labels = np.load("/home/rom/Dokumente/BA_Data/precomputed/10s_labels.npy")
access = np.load("/home/rom/Dokumente/BA_Data/precomputed/10s_access.npy")

print("Access:")
for a in access:
    print(a)

print(len(access))

for l in labels:
    print(l)

print(len(labels))
