#!/home/felix/miniconda3/bin/python
#$ -S /home/felix/miniconda3/bin/python
#$ -cwd
#$ -o /home/felix/digzyme/joblog/make_rxn_stdout/
#$ -e /home/felix/digzyme/joblog/make_rxn_stderr/


from datetime import datetime
from argparse import ArgumentParser


def append_mol(rxn_string, mol_dir, mol_list):
    for mol in mol_list:
        mol_path = f'{mol_dir}/{mol}.mol'
        with open(mol_path, encoding='utf8') as molfile:
            f_content = molfile.read()
        molfile.close()
        rxn_string += f'$MOL\n{f_content}'
    return rxn_string


def parser_setting():
    parser = ArgumentParser(
        prog='make_rxn.py',
        description=('Make RHEA-style RXN file'))
    parser.add_argument(
        'rheaid',
        action='store', type=str,
        help='RHEA reaction ID')
    parser.add_argument(
        '-m', '--moldir',
        action='store', type=str,
        help='Directory of MOL file(s)')
    parser.add_argument(
        '-o', '--outdir',
        action='store', type=str,
        help='Output file directory')
    parser.add_argument(
        '-r', '--reactants',
        action='store', type=str,
        help='List of reactant CHEBI(s)')
    parser.add_argument(
        '-p', '--products',
        action='store', type=str,
        help='List of product CHEBI(s)')
    parser.add_argument(
        '-d', '--delimiter',
        action='store', default=',',
        help='Delimiter of compound(s) list')
    args = parser.parse_args()
    return (parser, args)


def main(rheaid, moldir, outdir, reactants, products, delimiter):
    # Prepare variables
    timestamp = datetime.today().strftime('%d%m%Y%H%M')
    reactants = reactants.rstrip().split(delimiter)
    products = products.rstrip().split(delimiter)
    if len(reactants) > 0 and len(products) > 0:
        reactant_compounds = len(reactants)
        product_compounds = len(products)
    else:
        message = 'Please use \'-r\' and \'-p\' option'
        raise ValueError(message)

    # Write RXN content
    to_write = (f'$RXN\n\n'
                f'Felix  make_rxn{timestamp}  {rheaid}\n'
                f'RHEA:custom\n'
                f'  {reactant_compounds}  {product_compounds}\n')
    to_write = append_mol(to_write, moldir, reactants)
    to_write = append_mol(to_write, moldir, products)

    # Write RXN to output file
    output = f'{outdir}/{rheaid}.rxn'
    with open(output, mode='w', encoding='utf8') as outfile:
        outfile.write(to_write)
    outfile.close()


if __name__ == '__main__':
    parser, args = parser_setting()

    main(**vars(args))
