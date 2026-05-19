from Bio import SeqIO

record = SeqIO.read("memphis.gb", "genbank")

for feature in record.features:
    if feature.type == "mat_peptide":
        product = feature.qualifiers.get("product", [""])[0]
        seq = feature.extract(record.seq)
        print(f">{product}")
        print(seq.translate(to_stop=True))
