import gzip
import os
import sys


print("remove reads from paired fastq files that do not have excact given barcodes")
print("Easy but hard")
print("$1 = R1 fastq.gz file")
print("$2 = sample barcodes file")

barcodes = {}
with open(sys.argv[2]) as f:
    for line in f:
        words = line.split("\t")
        barcodes[words[0]] = words[1].rstrip()


if not os.path.exists("cleaned"):
    os.mkdir("cleaned")

r1File = sys.argv[1]
r2File = r1File.replace("_R1_", "_R2_")

r1 = gzip.open(r1File, "rt")
r2 = gzip.open(r2File, "rt")
r1Out = gzip.open("cleaned/" + r1File, "wt")
r2Out = gzip.open("cleaned/" + r2File, "wt")

sampleBarcode = barcodes[r1File]

print("sampleBarcode = " + sampleBarcode)

count = 0

while True:
    a1 = r1.readline()

    if not a1:
        break

    count += 1

    a2 = r1.readline()
    a3 = r1.readline()
    a4 = r1.readline()

    b1 = r2.readline()
    b2 = r2.readline()
    b3 = r2.readline()
    b4 = r2.readline()

    w1 = a1.split(" ")[1]
    barcode = w1.split(":")[3].rstrip()

    if barcode == sampleBarcode:
        r1Out.write(a1)
        r1Out.write(a2)
        r1Out.write(a3)
        r1Out.write(a4)

        r2Out.write(b1)
        r2Out.write(b2)
        r2Out.write(b3)
        r2Out.write(b4)
        

r1.close()
r2.close()
r1Out.close()
r2Out.close()


'''
for line in r1:
    count += 1

    if count > 10:
        break

    #if count % 4 == 1:
    #    print(line)

    print(line, end="")
'''






        