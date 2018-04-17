
def add_base(seq, pos, family, base):

    added = False

    while not added:
                            
        if base_pos in seq:

            if family_key in seq[base_pos]:

                if base in seq[base_pos][family_key]:
                    seq[base_pos][family_key][base] += 1

                else:
                    seq[base_pos][family_key][base] = 1

                added = True

            else:
                seq[base_pos][family_key] = {}

        else:
            seq[base_pos] = {}
            
