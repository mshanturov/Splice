import pandas as pd
import numpy

def make_df_column_str():
    column = set()
    stri = set()
    with open("C:/Users/liza_/Desktop/sww/other/distance.tsv", 'r') as f:
        line = f.readline()
        while line != '':
            if line[0] != "#":
                line_parts = line.split("\t")
                column.add(int(line_parts[1]))
                stri.add(int(line_parts[3]))
            line = f.readline()
        df = pd.DataFrame(index=stri, columns=column)
    return df, column, stri

def make_df_column_str2():
    df, column, stri = make_df_column_str()
    with open("C:/Users/liza_/Desktop/sww/other/distance.tsv", 'r') as r:
        line = r.readline()
        while line != '':
            if line[0] != "#":
                line_parts = line.split("\t")
                df[int(line_parts[1])][int(line_parts[3])] = float(line_parts[5])
            line = r.readline()
    df.to_excel("C:/Users/liza_/Desktop/sww/other/distance.xlsx")
make_df_column_str2()

# пока что херня
def make_martix_distance2():
    dict3 = {}
    dict5 = {}
    column = set()
    stri = set()
    with open("C:/Users/liza_/Desktop/sww/other/distances.tsv", 'r') as f:
        line = f.readline()
        while line != '':
            if line[0] != "#":
                line_parts = line.split("\t")
                column.add(int(line_parts[1]))
                stri.add(int(line_parts[3]))
            line = f.readline()
        column = sorted(column)
        stri = sorted(stri)
    count3 = 0
    count5 = 0
    for i in range(len(column)):
        if column[i] not in dict3:
            dict3[column[i]] = count3
            count3 += 1
    for j in range(len(stri)):
        if stri[j] not in dict5:
            dict5[stri[j]] = count5
            count5 += 1
    array = numpy.zeros((1443, 1443))
    with open("C:/Users/liza_/Desktop/sww/other/distances.tsv", 'r') as f:
        line = f.readline()
        while line != '':
            if line[0] != "#":
                line_parts = line.split("\t")
                array[dict3[int(line_parts[1])]][dict5[int(line_parts[3])]] = float(line_parts[4])
            line = f.readline()
    file = open("C:/Users/liza_/Desktop/sww/other/distances_test.txt", "w+")
    content = str(array)
    file.write(content)
    file.close()


