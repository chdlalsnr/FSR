import sys
from glob import glob


def main():

	flist = glob("*.csv")

	number = 0
	for fname in flist:

		f = open(fname,"read")
		lines = f.read().split("\n")

		number += len(lines)

	print(number)


if __name__ == "__main__":
        main()
