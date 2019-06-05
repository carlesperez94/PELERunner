def solve_ters_of_atom_section(pdb_lines_list):
    resid_saved = set()
    resid_jumped = set()
    index_to_ter = []
    for n, line in enumerate(pdb_lines_list):
        if line.startswith("ATOM"):
            resid = int(line[22:26])
            resid_saved.add(resid)
            res_previous = resid - 1
            if not res_previous in resid_saved and resid not in resid_jumped:
                resid_jumped.add(resid)
                index_to_ter.append(n)
    # Delete the first element that will be always detected as false positive
    index_to_ter.pop(0)
    # Reverse the order to insert things at the end of the list
    index_to_ter = sorted(index_to_ter, reverse=True)
    for index in index_to_ter:
        pdb_lines_list.insert(index, "TER\n")
    return pdb_lines_list


def solve_ters_and_water(pdb_lines_list):
    hoh_counter = 0
    h_counter = 1
    indexes_to_insert = []
    water_indexes = []
    for n, line in enumerate(pdb_lines_list):
        if line.startswith("HETATM"):
            resname = line[17:20].strip()
            if resname == "HOH":
                water_indexes.append(n)
                atom_name = line[12:16].strip()
                if "O" in atom_name:
                    atom_name = " OW "
                if "H" in atom_name:
                    hcheck = h_counter % 2
                    if hcheck == 1:
                        atom_name = "1HW "
                    elif hcheck == 0:
                        atom_name = "2HW "
                    h_counter += 1
                index = n - 1
                line = list(line)
                line[12:16] = atom_name
                line = "".join(line)
                pdb_lines_list[n] = line
                check_triad = hoh_counter % 3
                if check_triad == 1:
                    indexes_to_insert.append(index)
                hoh_counter += 1
    indexes_to_insert.append(water_indexes[-1]+1)
    indexes_to_insert = sorted(indexes_to_insert, reverse=True)
    for index in indexes_to_insert:
        if not pdb_lines_list[index-1].startswith("TER"):
            pdb_lines_list.insert(index, "TER\n")
    return pdb_lines_list


def add_ters_to_pdb(pdb_file):
    with open(pdb_file) as pdb:
        pdb_lines = pdb.readlines()
    atom_solved = solve_ters_of_atom_section(pdb_lines)
    water_solved = solve_ters_and_water(atom_solved)
    pdb_corrected = "".join(water_solved)
    return pdb_corrected

