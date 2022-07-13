def isochornost_cg(line):
    N = len(line.split("\n")[0])
    count_c, count_g = 0, 0
    for i in range(N):
        if line[i].lower().__eq__("c"):
            count_c += 1
        elif line[i].lower().__eq__("g"):
            count_g += 1
    return (count_c + count_g) / N


def isochornost_file():
    file1 = open("C:/Users/liza_/Desktop/sww/res3_upper.fasta", 'r')
    file2 = open("C:/Users/liza_/Desktop/sww/isochronost3.txt", 'w')
    lines = file1.readlines()
    for line in lines:
        if line[0] != ">":
            file2.write(str(isochornost_cg(line))+"\n")
isochornost_file()

