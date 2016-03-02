#!/usr/bin/env python
# author          : Pavel
# date            : 12.02.16
# version         : 0.1
# python_version  : 3
# copyright       : Pavel Polishchuk 2016
# license         : GPL3
#==============================================================================

import os
import argparse
from datetime import datetime
from itertools import combinations
from indigo import Indigo, IndigoException

indigo = Indigo()

"""
all frag_* functions return dict
general view: {fragment_name_1: (counter_fragment_1, counter_fragment_2), ...}
murcko returns only one fragment: {fragment_name_1: (counter_fragment_1)}
"""


def frag_rings(mol):
    output = dict()
    rings = []
    for ring in mol.iterateRings(1, mol.countAtoms()):
        atom_ids = [atom.index() for atom in ring.iterateAtoms()]
        smi = mol.createSubmolecule(atom_ids).canonicalSmiles().split(" ")[0]
        rings.append([smi, set(atom_ids)])
    rings = sorted(rings, key=lambda x: len(x[1]), reverse=True)
    for i, r in enumerate(rings):
        # if ring is not present in a bigger one than add counter-fragment to the output
        if all([not r[1].issubset(rings[j][1]) for j in range(i)]):
            if r[0] in output.keys():
                output[r[0]].append(mol.createSubmolecule(tuple(set(range(mol.countAtoms())).difference(r[1]))))
            else:
                output[r[0]] = [mol.createSubmolecule(tuple(set(range(mol.countAtoms())).difference(r[1])))]
    for k, v in output.items():
        output[k] = tuple(v)
    return output


def frag_murcko(mol):

    def add_atoms(atom, remove_atoms):
        for a in atom.iterateNeighbors():
            d = set([ai.index() for ai in a.iterateNeighbors()]).difference(remove_atoms)
            if len(d) == 1:
                remove_atoms.append(a.index())
                add_atoms(a, remove_atoms)

    # return Murcko framework
    remove_atoms = []
    for a in mol.iterateAtoms():
        if a.degree() == 1:
            remove_atoms.append(a.index())
            add_atoms(a, remove_atoms)

    # if all atoms were removed return empty dict
    if len(remove_atoms) == mol.countAtoms():
        return dict()
    else:
        return {mol.createSubmolecule(tuple(set(range(mol.countAtoms())).difference(remove_atoms))).canonicalSmiles().split()[0]: tuple([mol.createSubmolecule(remove_atoms)])}


def load_queries(query_fname):
    output = []
    for qmol in indigo.iterateSmilesFile(query_fname):
        m = indigo.loadQueryMolecule(qmol.smiles())
        if qmol.name():
            m.setName(qmol.name())
        else:
            m.setName(m.canonicalSmiles().split()[0])
        output.append(m)
    return output


def frag_query(mol, query_mols):
    output = dict()
    for qmol in query_mols:
        # several maps of a query mol can be
        res = []
        matcher = indigo.substructureMatcher(mol)
        for match in matcher.iterateMatches(qmol):
            ids = []
            for atom in qmol.iterateAtoms():
                if match.mapAtom(atom) is not None:
                    ids.append(match.mapAtom(atom).index())
            res.append(mol.createSubmolecule(tuple(set(range(mol.countAtoms())).difference(ids))))
        output[qmol.name()] = tuple(res)
    return output


def frag_auto(mol, query, max_cuts):

    def fix_cansmi_attach_point(smi):
        # fix cansmi name by replacement of attachment points: Cl%91.[*:1]%91 Cl[*:1]
        p = smi.split(".")
        output = p[0]
        for i in p[1:]:
            att = i.split('%')
            output = output.replace('%' + att[1], "(" + att[0] + ")")
        return output

    def frag_mol_by_cuts(mol, cut_list):

        output = []

        for cut in cut_list:
            a1 = mol.getAtom(cut[0])
            for nei in a1.iterateNeighbors():
                if nei.index() == cut[1]:
                    nei.bond().remove()
                    break
            # add attachment points
            mol.getAtom(cut[0]).setAttachmentPoint(1)
            mol.getAtom(cut[1]).setAttachmentPoint(1)

        for comp in mol.iterateComponents():
            cansmi = comp.clone().canonicalSmiles().split(" ")[0]
            cansmi = fix_cansmi_attach_point(cansmi)
            ids = [a.index() for a in comp.iterateAtoms()]
            output.append([cansmi, ids])

        return output

    def filter_dupl(frag_list):
        # input format: [['F', [0]], ['C#N', [3, 4]], ... ]
        res = dict()
        for item in frag_list:
            res[frozenset(sorted(item[1]))] = item[0]
        return [[v, list(k)] for k, v in res.items()]

    output = []

    matcher = indigo.substructureMatcher(mol)
    all_cuts = []
    for match in matcher.iterateMatches(query):
        ids = []
        for atom in query.iterateAtoms():
            if match.mapAtom(atom) is not None:
                ids.append(match.mapAtom(atom).index())
        all_cuts.append(ids)

    for i in range(1, max_cuts + 1):
        for comb in combinations(all_cuts, i):
            output.extend(frag_mol_by_cuts(mol.clone(), comb))

    # format: [['F', [0]], ['C#N', [3, 4]], ... ]
    output = filter_dupl(output)

    d = dict()
    for item in output:
        if item[0] not in d.keys():
            d[item[0]] = [mol.createSubmolecule(tuple(set(range(mol.countAtoms())).difference(item[1])))]
        else:
            d[item[0]].append(mol.createSubmolecule(tuple(set(range(mol.countAtoms())).difference(item[1]))))

    return d


def remove_single_h(frag_mols):
    output = dict()
    for frag_name, mols in frag_mols.items():
        curated_mols = []
        for mol in mols:
            hydrogens = [atom.index() for atom in mol.iterateAtoms() if atom.atomicNumber() == 1 and atom.degree() == 0]
            if len(hydrogens) < mol.countAtoms():
                curated_mols.append(mol.createSubmolecule(tuple(set(range(mol.countAtoms())).difference(hydrogens))))
        if curated_mols:
            output[frag_name] = tuple(curated_mols)
    return output


def remove_multicomponenet(frag_mols):
    output = dict()
    for frag_name, mols in frag_mols.items():
        output_mols = []
        for mol in mols:
            if mol.countComponents() == 1:
                output_mols.append(mol)
        if output_mols:
            output[frag_name] = tuple(output_mols)
    return output


def add_hydrogens(frag_mols):
    output = dict()
    for frag_name, mols in frag_mols.items():
        # output_mols = []
        for mol in mols:
            for atom in mol.iterateAtoms():
                for _ in range(atom.countImplicitHydrogens()):
                    a = mol.addAtom('H')
                    atom.addBond(a, 1)
        #     output_mols.append(mol)
        # if output_mols:
        #     output[frag_name] = tuple(output_mols)
    return frag_mols



def main_params(in_fname, title, out_fname, frag_scheme, terminal_only, remove_h, attach_h, keep, sep, query_fname, bond_smarts, max_cuts, verbose, error_fname):

    # remove output file if it exists
    if os.path.isfile(out_fname):
        os.remove(out_fname)

    input_format = os.path.splitext(in_fname)[1].lower()
    output_format = os.path.splitext(out_fname)[1].lower()

    if input_format == '.sdf':
        reader = indigo.iterateSDFile(in_fname)
    elif input_format == '.smi':
        reader = indigo.iterateSmilesFile(in_fname)
    else:
        print("Wrong input file extension.")
        exit()

    writer = indigo.writeFile(out_fname)

    if frag_scheme == 'query':
        query_mols = load_queries(query_fname)

    if frag_scheme == 'auto':
        bond_smarts = indigo.loadSmarts(bond_smarts)

    for i, mol in enumerate(reader):

        try:
            # get input mol name
            if input_format == '.sdf' and title is not None:
                mol_name = mol.getProperty(title)
            else:
                if mol.name():
                    mol_name = mol.name()
                else:
                    mol_name = 'MolID_' + str(i)

            if verbose:
                print(mol_name)

            if frag_scheme == 'rings':
                fragmented_mols = frag_rings(mol)
            elif frag_scheme == 'murcko':
                fragmented_mols = frag_murcko(mol)
            elif frag_scheme == 'query':
                fragmented_mols = frag_query(mol, query_mols)
            elif frag_scheme == 'auto':
                fragmented_mols = frag_auto(mol, bond_smarts, max_cuts)
            else:
                fragmented_mols = None

            # remove single detached hydrogens remaining after fragment removal
            if remove_h:
                fragmented_mols = remove_single_h(fragmented_mols)

            # remove all multi-component counter-fragments
            if terminal_only:
                fragmented_mols = remove_multicomponenet(fragmented_mols)

            # cap free valence with explicit hydrogenss
            if attach_h:
                fragmented_mols = add_hydrogens(fragmented_mols)

            # of no fragments were generated go to the next molecule
            if not fragmented_mols:
                continue

            # save input molecules if needed
            if keep:
                if output_format == '.sdf':
                    writer.sdfAppend(mol)
                elif output_format == '.smi':
                    open(out_fname, 'at').write(mol.smiles() + '\t' + mol_name + '\n')

            # save counter-fragments
            for frag_name, frag_mols in fragmented_mols.items():
                for frag_mol in frag_mols:
                    if frag_mol.countAtoms() == 0:
                        continue
                    # set names of output
                    name = mol_name + sep + frag_name
                    frag_mol.setName(name)
                    if output_format == '.sdf':
                        if title is not None:
                            frag_mol.setProperty(title, name)
                        writer.sdfAppend(frag_mol)
                    elif output_format == '.smi':
                        open(out_fname, 'at').write(frag_mol.smiles() + '\t' + name + '\n')

        except IndigoException as e:
            print('%s was skipped due to error' % mol.name())
            print(e)
            with open(error_fname, 'at') as f_err:
                f_err.write('%s\t%s\t%s\t%s\n' % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                                  os.path.abspath(in_fname), mol.name(), e))


def main():

    parser = argparse.ArgumentParser(description='Create structures with removed fragments of choice '
                                                 '(Indigo based implementation).')
    parser.add_argument('-i', '--in', metavar='input.sdf', required=True,
                        help='input sdf or SMILES file.')
    parser.add_argument('-t', '--title', metavar='sdf_title_field', default=None,
                        help='(optional) sdf field name containing titles of compounds, if it will omitted molecular '
                             'titles from sdf file will be used, if they are empty titles will be generated '
                             'automatically.')
    parser.add_argument('-o', '--out', metavar='output.sdf', required=True,
                        help='output sdf or smi (smiles) file.')
    parser.add_argument('-f', '--fragmentation', metavar='auto|rings|murcko|query', required=True,
                        help='fragmentation scheme - auto: rule-based, provide SMARTS to match cleaved bonds during '
                             'fragmentation; rings: all rings (fused rings considered as one ring system); '
                             'murcko: get murcko framework from each molecule; user: provide SMARTS patterns for '
                             'searched fragments.')
    parser.add_argument('-z', '--save_terminal_frags_only', action='store_true', default=False,
                        help='save only those fragments (terminal fragments) which produce just one counter-fragments. '
                             'Because scaffolds produce multi-component counter-fragments and that can lead to '
                             'problems for descriptors generation for some software. Default: false.')
    parser.add_argument('-a', '--attach_hydrogens', action='store_true', default=True,
                        help='automatically attach hydrogens to cap free valences after fragmentation. Default: true.')
    parser.add_argument('-x', '--remove_hydrogens', action='store_true', default=True,
                        help='remove single detached hydrogens remained after removal of a fragment. Default: true.')
    parser.add_argument('-k', '--keep_input_structures', action='store_true', default=True,
                        help='save input molecules in output file to simplify calculation of contributions.')
    parser.add_argument('-s', '--sep', metavar='separator_string', default='###',
                        help='separator which will be used for creation names of counter-fragments and will be '
                             'further used for splitting names (collisions with molecules and fragments names must '
                             'be avoided. Default : ###).')
    parser.add_argument('-q', '--query', metavar='user.smarts', default=None,
                        help='file containing fragments in SMARTS/SMILES format. '
                             'At each line the file contains tab-separated fields: 1) SMARTS/SMILES, 2) fragment name '
                             '(optional), 3) comma-separated list of atom numbers for recursive addition (optional). '
                             'Comment lines start with #')
    parser.add_argument('-r', '--rule', metavar='smarts_bond_pattern', default='[#6+0;!$(*=,#[!#6])]!@!=!#[*]',
                        help='SMARTS pattern for bond cleavage.')
    parser.add_argument('-u', '--max_cuts_number', metavar='integer', default=3,
                        help='maximal number of bonds cleaved simultaneously in rule-based fragmentation. Default: 3')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='show progress on the screen.')
    parser.add_argument('-e', '--error_file', metavar='log_file_name.txt', default="indigo_errors.txt",
                        help='save names of molecules which cause error to a text log file. Default file name '
                             'indigo_errors.txt in the directory of the input file.')

    args = vars(parser.parse_args())
    main_params(in_fname=args['in'],
                title=args['title'],
                out_fname=args['out'],
                frag_scheme=args['fragmentation'],
                attach_h=args['attach_hydrogens'],
                query_fname=args['query'],
                bond_smarts=args['rule'],
                max_cuts=args['max_cuts_number'],
                verbose=args['verbose'],
                error_fname=args['error_file'],
                keep=args['keep_input_structures'],
                sep=args['sep'],
                remove_h=args['remove_hydrogens'],
                terminal_only=args['save_terminal_frags_only'])


if __name__ == '__main__':
    main()
