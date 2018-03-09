import re

# open file and get data
f = open("./data_0216_9.txt")
data = f.read()
f.close()

# extract number from file data
numbers = re.findall('-?[0-9]+', data)

# string list -> int list
xy_list = list(map(int, numbers))



# delete serial number
s = 1
length = len(xy_list)

for i in reversed(range(length)):
    if i == 0:
        xy_list.pop(i)
    
    elif xy_list[i] == s and xy_list[i - 1] == -1:
        xy_list.pop(i)
        
        s -= 1


# delete area value
length = len(xy_list)

for i in reversed(range(length)):
    if i%3 == 2:
        xy_list.pop(i)




# substitute tuple (x, y) for list "pos"
x = 0
y = 0
pos = []

for i in range(len(xy_list)):
    if i%2 == 0:
        x = xy_list[i]
    else:
        y = xy_list[i]
        
        pos.append((x, y))


# file write
f = open("./_data_0216_9.txt", "w")

for i in range(len(pos)):
    x = pos[i][0]
    y = pos[i][1]
    
    if x != -1:
        f.write(str(x)+","+str(y))
        f.write("\n")
    else:
        f.write("\n")

f.close()









