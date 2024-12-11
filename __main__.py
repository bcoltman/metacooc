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

  MetaCoOc v%s

  Main dish:
    cooccurrence -> determine cooccurence based on detection of organisms in shotgun metagenomes

  Sides:
    setup -> download tables containing the sample metadata and presence/absence data
    filter -> filter the data by keyword or sample list
    ratio -> calculate (co)occurrence
    search -> query the sample metadata and presence/absence data tables
    list -> cook up a list of sample accessions based on filtering criteria

  Use: metacooc <command> -h for command specific help
    ''' % __version__)


def main():
    args = None
    if len(sys.argv) == 1:
        print_help()
        sys.exit(0)
    elif sys.argv[1] in {'-v', '--v', '-version', '--version'}:
        print(f"aastk: version {__version__} {__copyright__} {__author__}")
        sys.exit(0)
    elif sys.argv[1] in {'-h', '--h', '-help', '--help'}:
        print_help()
        sys.exit(0)
    else:
        parse_cli()


if __name__ == "__main__":
    main()
