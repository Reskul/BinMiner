if __name__ == '__main__':
    filepath_old = "C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\Bachelor_Arbeit\\illumina_10s_denovo\\illumina_10species.1.fq"
    filepath = "C:\\Users\\resku\\OneDrive\\Dokumente\\10SoSe21\\GCF_000001405.40_GRCh38.p14_genomic.fna"

    file = open(filepath, 'r')
    flag = True
    count = 0
    while flag:
        try:
            if file.readline().__contains__('>'):
                count += 1
                print("count++",count)
        except EOFError as e:
            flag = False
            print(e)
