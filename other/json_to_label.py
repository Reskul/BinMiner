import json
import matplotlib.pyplot as plt
from collections import Counter

if __name__ == '__main__':
    with open(f'f:\\10s.json', 'r') as json_file:
        label_dict: dict = json.load(json_file)
        keys = label_dict.keys()
        for key in keys:
            plt.clf()
            reads = label_dict[key]
            x = Counter(reads).keys()
            y_abs = list(Counter(reads).values())
            y_rel = [(v/sum(y_abs))*100 for v in y_abs]

            plt.title(key)
            plt.bar(x, y_rel)
            plt.xticks(rotation=60, ha='right')
            plt.tight_layout()
            plt.show()
            input()
