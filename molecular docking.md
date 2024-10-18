We used Autodock vina to perform the molecular docking.

```bash
vina_1.2.5_linux_x86_64 --ligand ligand.pdbqt --receptor receptor.pdbqt \
--scoring vina --spacing 1 --center_x x --center_y y --center_z z \
--size_x 126 --size_y 126 --size_z 126 --exhaustiveness 16 --seed 99 \
--out vinaResult.pdbqt --cpu 1
```

We used pymol to remove the heteroatoms in pdb file and located the center coordinate of the structure.
