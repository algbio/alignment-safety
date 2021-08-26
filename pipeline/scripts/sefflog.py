import os
import glob

def main():
    logs = sorted(glob.glob("logs/stats/*.out"))
    for log in logs:
        jid = log.split("/")[-1].split(".")[0]
        os.system(f"seff {jid} >> {log}")

if __name__ == '__main__':
    main()