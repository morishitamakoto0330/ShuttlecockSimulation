import re

# open file and get data
f = open('data4.txt')
data = f.read()
f.close()

# extract number from file data
numbers = re.findall('[0-9]+', data)

# delete area value
length = len(numbers)
for i in reversed(range(length)):
    if i%3 == 2:
        numbers.pop(i)

# string list -> int list
xy_list = list(map(int, numbers))

print(xy_list)

