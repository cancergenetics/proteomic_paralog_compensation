import argparse


def parse_args():
    """
    Parse command line arguments for running paralog compensation and collateral loss analysis.
    """
    parser = argparse.ArgumentParser(
        description='Run paralog compensation and collateral loss analysis.'
    )
    parser.add_argument(
        '--runwith',
        choices=['trans', 'prot', 'prot_residual'],
        required=True,
        help='Run pipeline with mRNA data, proteomic data, or protein residuals.'
    )
    parser.add_argument(
        '--nb_workers',
        type=int,
        required=False,
        default=1,
        help='Number of workers for parallel processing- no parallel if unspecified.'
    )
    parser.add_argument(
        '--HAP1_overlap',
        action="store_true",
        required=False,
        default=False,
        help='If set, restrict to pairs tested in HAP1 overlap analysis.'
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()