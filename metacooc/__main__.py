#!/usr/bin/env python3
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import sys
import traceback

from metacooc import __author__, __copyright__, __version__
from metacooc.cli import parse_cli

def print_help():
    print('''\

  metacooc v%s

  Starter:
    download -> Download required files with metagenome metadata and taxonomic profiles.

  Main dishes:
    cooccurrence -> Identify co-occurring prokaryotes across shotgun metagenomes.
    association  -> Identify prokaryotes associated with metadata terms, biomes, or other taxa.
    structure    -> Quantify community structure (e.g. nestedness, overlap, checkerboard patterns) in taxonomic profiles.

  Sides:
    search -> Query the metagenome metadata and taxonomic profiles.
    filter -> Filter the taxonomic profiles by keyword or sample list.
    analysis -> Perform association, cooccurrence or structural analysis.
    format -> Format sandpiper output, or other taxonomic profiles, to metacooc input.
    plot -> Generate plots based on the output of other menu items.
    biome_distribution -> Output a TSV file of taxa occurrence per pre-determined biomes.

  Use: metacooc <command> -h for command specific help
    ''' % __version__)


def main():
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f"metacooc: version {__version__} {__copyright__} {__author__}")
        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    elif sys.argv[1] not in {'download', 'cooccurrence', 'association', 'structure', 'search', 'filter', 'analysis', 'format', 'plot', 'biome_distribution'}:
        print(f"program not on the menu, choose from the options listed below ")
        print_help()
        sys.exit(0)
    else:
        parse_cli()


if __name__ == "__main__":
    main()
