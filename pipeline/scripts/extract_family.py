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
            scop_class = data[3]
            scop_fold = data[4]
            scop_superfamily = data[5]
            scop_family = data[6]
            nr_accession = data[1]
            family = scop_class + "-" + scop_fold + "-" + scop_superfamily + "-" + scop_family
            if not family in families:
                families[family] = []
            families[family].append(nr_accession)
    finally:
        f.close()
        m = 0
        m_i = ""
        # Largest family
        for family in families.keys():
            print(f"{family} : {len(families[family])}")
            if m < len(families[family]):
                m = len(families[family])
                m_i = family

        print(f"{m_i} : {m}")

    # Write family to file
    # for family in families.keys():
    # fam = "a-118-1-27"
    # with open("./data/scop_families/family_" + fam, "a") as file:
    #     for nr_accession in families[fam]:
    #         file.write(nr_accession + "\n")


def main():
    if(len(sys.argv) > 0):
        read_clusters(sys.argv[1])


if __name__ == '__main__':
    main()