'''
f = open('./_data_0216_9.txt', 'r')
lines = f.read().split('\n')
lines.pop(len(lines) - 1)
f.close()

xy_list = []

for line in lines:
    s = line.split(',')
    xy_list.append(list(map(int, s)))

v = 0
prev_y = 0
y = 0
sum_v = 0
count = 0


for index, xy in enumerate(xy_list):
    print(xy)

    y = xy[1]
    if index == 0:
        pass
    else:
        v = (y - prev_y)*60*4.08/1000;
        sum_v += v
        count += 1
        print("{0} [m/s]".format(v))

    prev_y = y

print("count : {0}".format(count))
print("average velocity : {0}".format(sum_v/count))
'''

print("terminal velocity={0}[m/s]".format(711.0*60.0*4.08/1000.0/22.0))


