import gemmi

inp = "3ZP0.cif"
out = "3ZP0_glycan_A.pdb"

target_chains = {"A"}

st = gemmi.read_structure(inp)

new = gemmi.Structure()
new.name = st.name
new.spacegroup_hm = st.spacegroup_hm
new.cell = st.cell

for model in st:
    new_model = gemmi.Model(model.num)

    for chain in model:
        new_chain = gemmi.Chain(chain.name)

        if chain.name in target_chains:
            for res in chain:
                new_chain.add_residue(res.clone())

        if len(new_chain) > 0:
            new_model.add_chain(new_chain)

    if len(new_model) > 0:
        new.add_model(new_model)

new.write_pdb(out)
print(f"Wrote {out}")

for model in new:
    for chain in model:
        print(f"Chain {chain.name}: {len(chain)} residues")
        for res in chain:
            print(f"  {res.name} {res.seqid}")
