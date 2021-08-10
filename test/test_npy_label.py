import numpy as np

labels = np.load("/home/rom/Dokumente/BA_Data/precomputed/10s_labels.npy")
access = np.load("/home/rom/Dokumente/BA_Data/precomputed/10s_access.npy")

for l in access:
    print(l)

print(len(access))
