#!/usr/bin/env python3
import sys
import re

GTF_FILE = "Hsalinarum_TSS.gtf"   # ajuste se necessário
WIN = 250
OUT = "promotor_gene.bed"

with open(GTF_FILE) as f, open(OUT, "w") as out:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        if parts[2] != "TSS":
            continue
        # Parseia atributos
        attrs = {}
        for attr in parts[8].split(";"):
            attr = attr.strip()
            if " " in attr:
                key, val = attr.split(" ", 1)
                val = val.strip('"')
                attrs[key] = val
            elif "=" in attr:
                key, val = attr.split("=", 1)
                attrs[key] = val.strip('"')
        tss_class = attrs.get("tss_class", "")
        if tss_class != "P":
            continue
        # Extrai gene de tss_location
        gene = ""
        loc = attrs.get("tss_location", "")
        if loc:
            m = re.search(r'gene_([^;]+)', loc)
            if m:
                gene = m.group(1)
            else:
                gene = "gene_" + str(start)
        if not gene:
            gene = "gene_" + str(start)
        # Janela 250 pb upstream
        if strand == "+":
            s = start - 1 - WIN
            e = start - 1
        else:
            s = start
            e = start + WIN
        if s < 0:
            s = 0
        out.write(f"{chrom}\t{s}\t{e}\t{gene}\t.\t{strand}\n")

print(f"Arquivo {OUT} gerado.")
