import argparse
import time


def get_time_ms():
    return round(time.time() * 1000)


def main(args):
    if args.action == "a":
        pass
    elif args.action == "i":
        pass
    elif args.action == "p":
        pass
    elif args.action == "w":
        pass
    else:
        pass
# Returns true if arguments are as required
def check_args(parser, args):
    if args.action == "a":
        return True
    elif args.action == "i":
        return True
    elif args.action == "p":
        parser.error("Action \"p\" requires --id")
        return "id" in vars(args)
    elif args.action == "w":
        parser.error("Action \"w\" requires --id")
        return "id" in vars(args)
    else:
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("action", type=str, help="a/i/p/w")
    parser.add_argument("db", type=str, help="database in fasta-format")
    parser.add_argument("clusters", type=str, help="cluster file")
    parser.add_argument("--min", type=int, default=10, help="minimum size of the cluster {10}")
    parser.add_argument("--max", type=int, default=1000, help="maximum size of the cluster {1000}")
    parser.add_argument("--id", type=str, help="cluster id to separate")
    args = parser.parse_args()
    if check_args(parser, args):
        main(args)
