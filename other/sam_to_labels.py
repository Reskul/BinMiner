import numpy as np
import json
from collections import Counter

if __name__ == '__main__':
    flag = False
    cnt = 0
    contigs = np.empty(0)
    labels_dict = {}
    filename = '10s'
    THRESHOLD = 10.0  # Threshold in percent to assign an organism to a contig
    # Test
    # with open('f:\\10s.sam', 'r') as sam_:
    #     test_file = open('f:\\10s_part.sam', 'w')
    #     for i in range(5000):
    #         l = sam_.readline()
    #         test_file.write(l)
    with open(f'f:\\{filename}.sam', 'r') as sam:
        for line in sam:
            if flag:
                content = line.split('\t')
                # print(content[0].split(':')[0], content[2])
                if content[2] not in labels_dict:
                    labels_dict[content[2]] = [content[0].split('.')[0]]
                else:
                    labels_dict[content[2]].append(content[0].split('.')[0])
            if line.startswith('@PG'):
                flag = True
                print(f'{cnt} contigs counted.')
            if line.startswith('@SQ'):
                contigs = np.append(contigs, line.split(':')[1].split('\t')[0])
                cnt += 1
        contigs = contigs[1:]
        print("Saving contigs list...")
        with open(f'f:\\10s_contigs_list.txt', 'w') as contigs_out:
            contigs_out.writelines(contigs)
        print("Saving to json...")
        export_file = open(f"f:\\{filename}.json", "w")
        export_file.write(json.dumps(labels_dict))
        export_file.close()

        keys = list(labels_dict.keys())
        idx = keys.index('*')
        keys.pop(idx)
        print("Permutation? ", Counter(keys) == Counter(contigs))

        # labels = np.empty((len(contigs), 2))
        assignments = []
        for key in keys:
            reads = labels_dict[key]
            x = list(Counter(reads).keys())
            y_abs = list(Counter(reads).values())
            y_rel = [(v / sum(y_abs)) * 100 for v in y_abs]
            max_val = max(y_abs)
            max_idx = y_abs.index(max_val)
            print(y_rel[max_idx])
            if y_rel[max_idx] >= THRESHOLD:
                assignments.append(x[max_idx])
            else:
                assignments.append('TNM_ERROR')  # Threshold Not Met Error
        print("Saving labels...")
        with open(f'f:\\10s_labels.txt', 'w') as out:
            i = 0
            for key in keys:
                out.write(f"{key}\t{assignments[i]}\n")
                i += 1
