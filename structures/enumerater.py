def quat_angle(c, i, j, k, l):
    print i, c, j, 1
    print i, c, k, 1
    print i, c, l, 1
    print j, c, k, 1
    print j, c, l, 1
    print l, c, k, 1

def dihed_set(iset, j, k, lset):
    for i in iset:
        for l in lset:
            print i, j, k, l, 3
def pair(iset, jset):
    for i in iset:
        for j in jset:
            print i, j

quat_angle(1, 2, 3, 4, 5)
quat_angle(2, 1, 6, 7, 8)
quat_angle(3, 1, 9, 10, 11)
quat_angle(4, 1, 12, 13, 14)
quat_angle(5, 1, 15, 16, 17)

dihed_set([6, 7, 8], 2, 1, [3, 4, 5])
dihed_set([9, 10, 11], 3, 1, [2, 4, 5])
dihed_set([12, 13, 14], 4, 1, [2, 3, 5])
dihed_set([15, 16, 17], 5, 1, [2, 3, 4])




pair([6, 7, 8], [3, 4, 5])
pair([9, 10, 11], [2, 4, 5])
pair([12, 13, 14], [2, 3, 5])
pair([15, 16, 17], [2, 3, 4])




