if __name__ == '__main__':
    in_filename = "f:\\10s_coverage.txt"
    out_filename = "f:\\10s_depth.txt"
    in_file = open(in_filename, 'r')
    out_file = open(out_filename, 'w')

    print(in_file.readline())
    idx = input("Indices(',' delimited):").split(',')
    for line in in_file:
        content = line.split('\t')
        out_file.write(f"{content[int(idx[0].strip())]}\t{content[int(idx[1].strip())]}\n")
        out_file.flush()

    out_file.close()
    in_file.close()

