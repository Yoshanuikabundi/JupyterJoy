def make_martini_ndx(top, f=None):
    """Construct a gromacs ndx file from a martini MDTraj topology"""
    ndx = {}
    for res in top.residues:
        name = res.name.upper()
        indices = [a.index+1 for a in res.atoms]
        try:
            ndx[name] += indices
        except KeyError:
            ndx[name] = indices
    ndx['system'] = [a.index+1 for a in top.atoms]
    ndx['solvent'] = sorted(ndx.get('PW', []) 
                            + ndx.get('W', []) 
                            + ndx.get('WF', []) 
                            + ndx.get('NA+', []) 
                            + ndx.get('CL-', []) 
                            + ndx.get('NC3+', []) 
                            + ndx.get('CA+', [])
                           )
    ndx['non-solvent'] = [a.index+1 for a in top.atoms if a.index+1 not in ndx['solvent']]
    if f: 
        for k,v in ndx.items():
            f.write("[ {} ]\n".format(k))
            for n,i in enumerate(v):
                f.write("{: >6} ".format(i))
                if not (n+1)%5: f.write('\n')
            f.write("\n\n")
    return ndx