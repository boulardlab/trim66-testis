d = {}
with open("GRCm38.primary_assembly.genome.fa") as fin:
    for line in fin:
        if line.startswith(">"):
            k = line.strip()
            d[k] = ""
        else:
            d[k] += line.strip()

d = {k: v for k in d.items() if "_" not in k or "." not in k}

with open("GRCm38.primary_assembly.genome.filter.fa") as fout:
    for k in d:
        fout.write(f"{k}\n")
        fout.write("\n".join([d[k][i:i+60] for i in range(0, len(d[k]), 60)]))

