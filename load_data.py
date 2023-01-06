import re

def SRD101():
    file=open("./bond_data/Standard Reference Database 101.txt")
    lines=file.readlines()
    print(lines)

if __name__=="__main__":
    SRD101()