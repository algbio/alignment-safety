import sys

families = {}

def read_clusters(filename):
    f = open(filename, "r")
    try:
        while f:
            line = f.readline()
            if not line:
                break
            data = line.rstrip("\n").split("\t")
            family = data[6]
            nr_accession = data[1]
            if not family in families:
                families[family] = []
            families[family].append(nr_accession)
    finally:
        f.close()
        m = 0
        m_i = ""
        # Largest family
        for family in families.keys():
            if m < len(families[family]):
                m = len(families[family])
                m_i = family

        print(f"{m_i} : {m}")

    # Write families to file
    # for family in families.keys():
    #     with open("./data/scop_families/family_" + family, "a") as file:
    #         for nr_accession in families[family]:
    #             file.write(nr_accession + "\n")


def main():
    if(len(sys.argv) > 0):
        read_clusters(sys.argv[1])


if __name__ == '__main__':
    main()