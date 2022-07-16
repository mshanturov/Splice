import pandas as pd


def make_martix_distance():
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
        column = sorted(column)
        stri = sorted(stri)
        df = pd.DataFrame(index=stri, columns=column)
    with open("C:/Users/liza_/Desktop/sww/other/distance.tsv", 'r') as r:
        line = r.readline()
        while line != '':
            if line[0] != "#":
                line_parts = line.split("\t")
                df[int(line_parts[1])][int(line_parts[3])] = float(line_parts[5])
            line = r.readline()
    df.to_excel("C:/Users/liza_/Desktop/sww/other/distance.xlsx")
    print(df)

make_martix_distance()

# пока что херня
def make_martix_distance2():
    dict3 = {}
    dict5 = {}
    count3 = 0
    count5 = 0
    array = [[0]*1443]*1443
    with open("C:/Users/liza_/Desktop/sww/other/distances.tsv") as f:
        line = f.readline()
        while line != '':
            if line[0] != "#":
                line_parts = line.split("\t")
                if line_parts[1] not in dict3:
                    dict3[(line_parts[1])] = count3
                    count3 += 1
                if line_parts[3] not in dict5:
                    dict5[(line_parts[3])] = count5
                    count5 += 1
                array[int(dict3[(line_parts[1])])][int(dict5[(line_parts[3])])] = (line_parts[4])
            line = f.readline()
    print(array)


