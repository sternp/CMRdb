#!/usr/bin/env python
###############################################################################
# cmrdb.py - Info about cmrdb.py
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
from cmrdb.__init__ import __version__
import cmrdb.config.config as Config
__author__ = "Peter Sternes"
__copyright__ = "Copyright 2023"
__credits__ = ["Peter Sternes"]
__license__ = "GPL3"
__maintainer__ = "Peter Sternes"
__email__ = "peter.sternes near qut.edu.au"
__status__ = "Development"


#%%################################ System imports ###########################
import sys
import argparse
import logging
import os
from datetime import datetime
import subprocess

# Local imports
from snakemake import utils
from snakemake.io import load_configfile
from ruamel.yaml import YAML  # used for yaml reading with comments

# Debug
debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

#%%############################  Exceptions  ################################
#############################################################################

class BadTreeFileException(Exception):
    pass

################################  Functions  ################################
#############################################################################
def centerify(text, width=-1):
  lines = text.split('\n')
  width = max(map(len, lines)) if width == -1 else width
  return '\n'.join(line.center(width) for line in lines)


def phelp():
    print(
"""
           _____ __  __ _____     _ _     
          / ____|  \/  |  __ \   | | |    
         | |    | \  / | |__) |__| | |__  
         | |    | |\/| |  _  // _` | '_ \ 
         | |____| |  | | | \ \ (_| | |_) |
          \_____|_|  |_|_|  \_\__,_|_.__/  Centre for Microbiome Research, QUT
                                          
                                          
        A pipline for the assembly of biobank samples and merging into the CMR genome database.

        process - Full pipeline. Raw reads -> assembled, QCed, annotated genomes. Can skip the read QC.
        merge - [NOT BUILT YET] Merge a genome into the underlying CMRdb

        Type 'cmrdb {process,merge} --help' for specific information

"""
)


def main():

    ############################  Main Parser  ##############################
    #########################################################################
    main_parser = argparse.ArgumentParser(
        prog='cmrdb',
        formatter_class=CustomHelpFormatter,
        add_help=False)

    main_parser.add_argument(
        '--version',
        action='version',
        version=__version__,
        help='Show version information.'
        )

    main_parser.add_argument(
        '--verbosity',
        help='1 = critical, 2 = error, 3 = warning, 4 = info, 5 = debug. Default = 4 (logging)',
        type=int,
        default=4
        )

    main_parser.add_argument(
        '--log',
        help='Output logging information to file',
        default=False
        )

    main_parser.add_argument(
        '--conda-prefix',
        help='Path to the location of installed conda environments, or where to install new environments',
        dest='conda_prefix',
        default=Config.get_conda_path(),
        metavar='<path>'
    )

    main_parser.add_argument(
        '-n', '--n-cores',
        help='Maximum number of cores available for use.',
        dest='n_cores',
        default=8,
        metavar='<num>'
    )

    main_parser.add_argument(
        '-m', '--max-memory',
        help='Maximum memory for available usage in gigabytes, GB',
        dest='max_memory',
        default=64,
        metavar='<num>'
    )

    main_parser.add_argument(
        '-o', '--output',
        help='Output directory',
        dest='output',
        default='./',
        metavar='<dir>'
    )

    main_parser.add_argument(
        '--dry-run',
        help='Perform snakemake dry run, tests workflow order and conda environments',
        type=str2bool,
        nargs='?',
        const=True,
        dest='dryrun',
        default=False
    )

    main_parser.add_argument(
        '--conda-frontend',
        help='Which conda frontend to use',
        dest='conda_frontend',
        nargs=1,
        default="mamba",
        choices=["conda", "mamba"],
        metavar='<type>'
    )




    subparsers = main_parser.add_subparsers()


    ##########################  sub-parser  ###########################
    parser_process = subparsers.add_parser('process',
                                        parents=[main_parser],
                                        formatter_class=CustomHelpFormatter,
                                        add_help=True, description='''

                                ~ FULL PIPELINE - QC, assembly, annotation ~
    How to use process:

    cmrdb process
        -1 reads_R1.fastq \\
        -2 reads_R2.fastq \\
        --long nanopore.fastq \\    
        -n 24 \\
        -m 128 \\
        -a isolate \\    
        -o output_directory
    ''')

    parser_process.add_argument(
        '-1',
        help='Forward short FASTQ reads',
        dest='pe1',
        default="none",
        metavar='<file1>'
    )

    parser_process.add_argument(
        '-2',
        help='Reverse short FASTQ reads',
        dest='pe2',
        default="none",
        metavar='<file2>'
    )

    parser_process.add_argument(
        '--long',
        help='Long reads (Optional; Nanopore only)',
        dest='long',
        default="none",
        metavar='<file3>'
    )

    parser_process.add_argument(
        '--sequencer-source',
        help='Used for QC/trimming of short reads: NexteraPE, TruSeq2, TruSeq3, none.',
        dest='sequencer_source',
        default="TruSeq3",
        metavar='<type>'
    )

    parser_process.add_argument(
        '--skip-qc',
        help='Skip the read QC step',
        dest='skip_qc',
        action='store_true',
        default=False
    )

    parser_process.add_argument(
        '-w', '--workflow',
        help=argparse.SUPPRESS,
        dest='workflow',
        default='process'
    )
    
    parser_process.add_argument(
        '-a', '--assembly-mode',
        help='Adjust the SPAdes assembler according to sample type. Typical uses are metagenomic (meta), isolate (isolate) or single cell (sc) mode',
        dest='assembly_mode',
        choices=['isolate', 'meta', 'sc','bio','rna','plasmid'],
        default='isolate'
    )
    
    ############################## Parsing input ##############################
    if (len(sys.argv) == 1 or len(sys.argv) == 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        phelp()
    else:
        args = main_parser.parse_args()
        time = datetime.now().strftime('%H:%M:%S %d-%m-%Y')

        if args.log:
            if os.path.isfile(args.log):
                raise Exception("File %s exists" % args.log)
            logging.basicConfig(filename=args.log,
                                level=debug[args.verbosity],
                                format='%(asctime)s %(levelname)s: %(message)s',
                                datefmt='%m/%d/%Y %I:%M:%S %p')
        else:
            logging.basicConfig(level=debug[args.verbosity],
                                format='%(asctime)s %(levelname)s: %(message)s',
                                datefmt='%m/%d/%Y %I:%M:%S %p')
        logging.info("Time - %s" % (time))
        logging.info("Command - %s" % ' '.join(sys.argv))
        logging.info("Version - %s" % __version__)

        prefix = args.output
        if not os.path.exists(prefix):
            os.makedirs(prefix)

        #fill-in Namespace for attributes which only appear in specific subparsers
        params=['pe1', 'pe2', 'long', 'n_cores', 'max_memory', 'output', 'conda_prefix', 'sequencer_source','skip_qc','workflow', 'assembly_mode' ]
        for i in params:
            try:
                getattr(args, i)
            except AttributeError:
                e = 'args.' + i + ' = \'none\''
                exec(e)
            else:
                pass

        processor = cmrdb(args.pe1,
                                args.pe2,
                                args.long,
                                int(args.n_cores),
                                int(args.max_memory),
                                args.output,
                                args.conda_prefix,
                                args.sequencer_source,                      
                                args.skip_qc,
                                args.workflow,
                                args.assembly_mode,
                                args)

        processor.make_config()
        processor.run_workflow(workflow=args.workflow, cores=int(args.n_cores), dryrun=args.dryrun, conda_frontend=args.conda_frontend)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


def update_config(config):
    """
    Populates config file with default config values.
    And made changes if necessary.
    """

    # get default values and update them with values specified in config file
    default_config = make_default_config()
    utils.update_config(default_config, config)

    return default_config


#################################  Classes  ###################################
###############################################################################

class CustomHelpFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        return text.splitlines()

    def _get_help_string(self, action):
        h = action.help
        if '%(default)' not in action.help:
            if action.default != '' and \
               action.default != [] and \
               action.default != None \
               and action.default != False:
                if action.default is not argparse.SUPPRESS:
                    defaulting_nargs = [argparse.OPTIONAL,
                                        argparse.ZERO_OR_MORE]

                    if action.option_strings or action.nargs in defaulting_nargs:

                        if '\n' in h:
                            lines = h.splitlines()
                            lines[0] += ' (default: %(default)s)'
                            h = '\n'.join(lines)
                        else:
                            h += ' (default: %(default)s)'
        return h

    def _fill_text(self, text, width, indent):
        return ''.join([indent + line for line in text.splitlines(True)])

class cmrdb:
    def __init__(self,
                 pe1="none",
                 pe2="none",
                 long="none",
                 n_cores=8,
                 max_memory=32,
                 output=".",
                 conda_prefix=Config.get_conda_path(),
                 sequencer_source = "TruSeq3",
                 skip_qc=False,
                 workflow="none",
                 assembly_mode="isolate",
                 args=None
                 ):
        self.pe1 = pe1
        self.pe2 = pe2
        self.long = long
        self.threads = n_cores
        self.max_memory = max_memory
        self.output = output
        self.conda_prefix = conda_prefix
        self.sequencer_source = sequencer_source
        self.skip_qc = skip_qc
        self.workflow = workflow
        self.assembly_mode = assembly_mode

    def make_config(self):
        """
        Reads template config file with comments from ./config.yaml
        updates it by the parameters provided.
        Args:
            ADD SOME INFO HERE
        """

        self.config = os.path.join(self.output, 'config.yaml')

        yaml = YAML()
        yaml.version = (1, 1)
        yaml.default_flow_style = False

        template_conf_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                          "config.yaml")

        with open(template_conf_file) as template_config:
            conf = yaml.load(template_config)

        if self.pe1 != "none":
            self.pe1 = os.path.abspath(self.pe1)
        if self.pe2 != "none":
            self.pe2 = os.path.abspath(self.pe2)
        if self.long != "none":
            self.long = os.path.abspath(self.long)    
        if self.output != "none":
            self.output = os.path.abspath(self.output)
        if self.sequencer_source != "TruSeq3":
            self.sequencer_source = self.sequencer_source
        if self.skip_qc != False:
            self.skip_qc = self.skip_qc
        if self.workflow != "none":
            self.workflow = self.workflow
        if self.assembly_mode != "isolate":
            self.assembly_mode = self.assembly_mode    

        conf["short_reads_1"] = self.pe1
        conf["short_reads_2"] = self.pe2
        conf["long_reads"] = self.long
        conf["n_cores"] = self.threads
        conf["max_memory"] = self.max_memory
        conf["output"] = self.output
        conf["sequencer_source"] = self.sequencer_source
        conf["skip_qc"] = self.skip_qc
        conf["workflow"] = self.workflow
        conf["assembly_mode"] = self.assembly_mode


        with open(self.config, "w") as f:
            yaml.dump(conf, f)
        logging.info("Configuration file written to %s\n" % self.config)

    def validate_config(self):
        load_configfile(self.config)


    def run_workflow(self, workflow="process", cores=8, profile=None, dryrun=False, conda_frontend="mamba", snakemake_args = ""):
        """Runs the CMRdb pipeline
        By default all steps are executed
        Needs a config-file which is generated by given inputs.
        Most snakemake arguments can be appended to the command for more info see 'snakemake --help'
        """

        if not os.path.exists(self.config):
            logging.critical(f"config-file not found: {self.config}\n")
            sys.exit(1)

        self.validate_config()

#        conf = load_configfile(self.config)

        cmd = (
            "snakemake --snakefile {snakefile} --directory {working_dir} "
            "{jobs} --rerun-incomplete "
            "--configfile '{config_file}' --nolock "
            " {profile} {conda_frontend} --use-conda {conda_prefix} {dryrun} "
#            " {target_rule} "
            " {args} "
        ).format(
            snakefile=get_snakefile(),
            working_dir=self.output,
            jobs="--jobs {}".format(cores) if cores is not None else "",
            config_file=self.config,
            profile="" if (profile is None) else "--profile {}".format(profile),
            dryrun="--dryrun" if dryrun else "",
            args=" ".join(snakemake_args),
#            target_rule=workflow if workflow != "None" else "",
            conda_prefix="--conda-prefix " + self.conda_prefix,
            conda_frontend="--conda-frontend " + conda_frontend
        )
        logging.info("Executing: %s" % cmd)
        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError as e:
            # removes the traceback
            logging.critical(e)
            exit(1)

if __name__ == '__main__':

    sys.exit(main())
