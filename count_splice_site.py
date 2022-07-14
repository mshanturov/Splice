import csv

filename3 = "res3"
filename5 = "res5"
file3 = open("C:/Users/liza_/Desktop/sww/{}.fa".format(filename3), 'r')
file5 = open("C:/Users/liza_/Desktop/sww/{}.fa".format(filename5), 'r')
dict_splice3, dict_splice5 = {}, {}
count_splice3, count_splice5 = 0, 0
N = 5

lines3 = file3.readlines()
for line3 in lines3:
    if line3[0] != ">":
        count_splice3 += 1
        if line3[:N].upper() in dict_splice3:
            dict_splice3[line3[:N].upper()] = dict_splice3[line3[:N].upper()] + 1
        else:
            dict_splice3[line3[:N].upper()] = 1

lines5 = file5.readlines()
for line5 in lines5:
    if line5[0] != ">":
        count_splice5 += 1
        if line5[:N].upper() in dict_splice5:
            dict_splice5[line5[:N].upper()] = dict_splice5[line5[:N].upper()] + 1
        else:
            dict_splice5[line5[:N].upper()] = 1


dict_splice3 = dict(sorted(dict_splice3.items(), key=lambda item: item[1], reverse=True))
dict_splice5 = dict(sorted(dict_splice5.items(), key=lambda item: item[1], reverse=True))

with open('C:/Users/liza_/Desktop/sww/sequences_{}.csv'.format(filename3), 'w') as f:
    w = csv.DictWriter(f, dict_splice3.keys())
    w.writeheader()
    w.writerow(dict_splice3)

with open('C:/Users/liza_/Desktop/sww/sequences_{}.csv'.format(filename5), 'w') as f:
    w = csv.DictWriter(f, dict_splice5.keys())
    w.writeheader()
    w.writerow(dict_splice5)

with open('C:/Users/liza_/Desktop/sww/count_{}_and_{}.csv'.format(filename3, filename5), 'w') as f:
    w = csv.writer(f)
    w.writerow(['all', 'splice3', 'splice5'])
    w.writerow([count_splice3+count_splice5, count_splice3, count_splice5])
