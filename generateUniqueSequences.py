#! /usr/bin/python3

import random
import sys
alphabet = ['A', 'C', 'G', 'U']

def generateUniqueFastas(uniques, min, max, step = 100):
    result = []
    for x in range(uniques):
        for y in range(min, max+1, step):
            res = ""
            print(">" + str(y) + "_randomseq")
            for z in range(y):
                # if ((z % 150) == 0 and z!=0):
                #     res += "\n"
                res += alphabet[random.randint(0, 3)]
            print(res)

generateUniqueFastas(10, 200, 900, 100)
generateUniqueFastas(10, 1000, 6000, 1000)