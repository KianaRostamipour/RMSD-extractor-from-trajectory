"""
This code reads RMSD values from a trajectory analysis file and generates a bash script (generate_structures.sh). The generated bash script is used to extract PDB structures from trajectory files at specific time points corresponding to RMSD quartiles.
"""
# Replace 'PDB ID' with the appropriate protein name as needed. The PDB ID here is PDBID.
import numpy as np
#change the directories accordingly
rmsd_file_path = 'analysis/rmsd_all_segments.xvg'
bash_script_path = 'generate_structures.sh'
def read_rmsd_data(file_path):
    times, values = [], []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith(('@', '#')):
                parts = line.split()
                times.append(float(parts[0]))
                values.append(float(parts[1]))
    return np.array(times), np.array(values)
times, rmsd_values = read_rmsd_data(rmsd_file_path)
quartiles_indices = np.linspace(0, len(rmsd_values)-1, 9, endpoint=True).astype(int)
with open(bash_script_path, 'w') as bash_script:
    bash_script.write("#!/bin/bash\n\n")
    bash_script.write("base_dir=$(basename \"$(pwd)\")\n")
    bash_script.write("mkdir -p RMSD_quartile_pdb_${base_dir}\n")
    bash_script.write("for file in molecule_*.itp; do\n")
    bash_script.write("  ligand_name=$(basename \"$file\" | grep -Eo '[0-9]+' | head -1)\n")
    bash_script.write("  echo \"Processing ligand $ligand_name\"\n")    
    for i in range(8):
        time_min = times[quartiles_indices[i]:quartiles_indices[i+1]][rmsd_values[quartiles_indices[i]:quartiles_indices[i+1]].argmin()]
        time_max = times[quartiles_indices[i]:quartiles_indices[i+1]][rmsd_values[quartiles_indices[i]:quartiles_indices[i+1]].argmax()]
        
        bash_script.write(f'  echo -e "20" | gmx trjconv -f production/md1_PDBID_${{ligand_name}}.xtc -s tpr/md1_PDBID_${{ligand_name}}.tpr -o RMSD_quartile_pdb_${{base_dir}}/min_time_quartile_{i+1}_${{ligand_name}}.pdb -n ndx/md1_PDBID_${{ligand_name}}.ndx -dump {time_min * 1000}\n')
        bash_script.write(f'  echo -e "20" | gmx trjconv -f production/md1_PDBID_${{ligand_name}}.xtc -s tpr/md1_PDBID_${{ligand_name}}.tpr -o RMSD_quartile_pdb_${{base_dir}}/max_time_quartile_{i+1}_${{ligand_name}}.pdb -n ndx/md1_PDBID_${{ligand_name}}.ndx -dump {time_max * 1000}\n\n')
    
    bash_script.write("done\n")

print(f"done!")

