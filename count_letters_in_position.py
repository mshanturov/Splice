def count_letters_in_position():
    file1 = open("C:/Users/liza_/Desktop/sww/res3_44_upper.fa", 'r')
    file2 = open("C:/Users/liza_/Desktop/sww/result_position.txt", 'w')
    lines = file1.readlines()
    n = len(lines[1].split("\n")[0])
    arr_a, arr_g, arr_t, arr_c = [0] * n, [0] * n, [0] * n, [0] * n
    for line in lines:
        if line[0] != ">":
            for i in range(len(line)):
                if line[i].lower().__eq__("a"):
                    arr_a[i] += 1
                elif line[i].lower().__eq__("c"):
                    arr_c[i] += 1
                elif line[i].lower().__eq__("g"):
                    arr_g[i] += 1
                elif line[i].lower().__eq__("t"):
                    arr_t[i] += 1
    for i in range(-45, 8):
        file2.write(str(i)+"\t")
    file2.write("\n")
    file2.write("A\t" + str(arr_a))
    file2.write("\n")
    file2.write("C\t" + str(arr_c))
    file2.write("\n")
    file2.write("G\t" + str(arr_g))
    file2.write("\n")
    file2.write("T\t"+str(arr_t))

count_A_in_position()