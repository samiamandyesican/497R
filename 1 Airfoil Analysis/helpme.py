import sys

def readlines(infile_name: str) -> list[str]:
    with open(infile_name) as file:
        return file.readlines()


def main(infile_name: str):
    lines = readlines(infile_name)
    x_cor = []
    y_cor = []
    for line in lines:
        x, y = line.split()
        print(f'{x} \t {y}')


if __name__ == "__main__":
    main(sys.argv[1])