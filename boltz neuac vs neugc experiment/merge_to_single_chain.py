import gemmi

inp = "9W44_glycan_DGJ.pdb"
out = "9W44_glycan_single.pdb"

st = gemmi.read_structure(inp)

new = gemmi.Structure()
new.name = st.name
new.spacegroup_hm = st.spacegroup_hm
new.cell = st.cell

for model in st:
    new_model = gemmi.Model(model.num)
    new_chain = gemmi.Chain("A")

    res_num = 1
    for chain in model:
        for res in chain:
            cloned = res.clone()
            cloned.seqid.num = res_num
            cloned.seqid.icode = ' '
            new_chain.add_residue(cloned)
            res_num += 1

    if len(new_chain) > 0:
        new_model.add_chain(new_chain)

    if len(new_model) > 0:
        new.add_model(new_model)

new.write_pdb(out)
print(f"Wrote {out}")

print(f"Models: {len(new)}")
for model in new:
    print(f"  Model {model.num}: {len(model)} chains")
    for chain in model:
        print(f"    Chain {chain.name}: {len(chain)} residues")
        for res in chain:
            print(f"      {res.name} {res.seqid}")
