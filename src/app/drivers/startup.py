import os
from argparse import ArgumentParser, RawTextHelpFormatter

VERSION_NUMBER = '2.10'


def parser_file_exists(currentparser, fn):
    if not os.path.exists(fn):
        currentparser.error('Input file %s not found!' % fn)
    else:
        return fn


def populate_parser(drivers):
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--version', action='version', version=VERSION_NUMBER)
    subparsers = parser.add_subparsers(dest='subcommand')
    subparsers.required = True
    subparsermap = {}
    for driver in drivers:
        driver.set_options()
        cmd = driver.get_commandname()
        subparsermap[cmd] = subparsers.add_parser(
            cmd, conflict_handler='resolve', help=driver.get_commandhelp())
        subparsermap[cmd].set_defaults(func=driver.start)
        subparsermap[cmd].add_argument('--version', action='version',
                                       version=VERSION_NUMBER)
        options = [{k: v for k, v in opt.items()}
                   for opt in driver.get_options()]
        for argoptions in options:
            argoptions['dest'] = argoptions.pop('driverattr')
            clarg = argoptions.pop('clarg')
            if 'type' in argoptions and argoptions['type'] == 'file':
                argoptions['type'] = lambda x: parser_file_exists(parser, x)
            elif 'type' in argoptions and argoptions['type'] == 'pick':
                argoptions.pop('picks')
                argoptions.pop('type')
            if not 'required' in argoptions:
                argoptions['required'] = True
            if type(clarg) == list:
                subparsermap[cmd].add_argument(*clarg, **argoptions)
            else:
                subparsermap[cmd].add_argument(clarg, **argoptions)
    return parser


def start_msstitch(exec_drivers, sysargs):
    """Passed all drivers of executable, checks which command is passed to
    the executable and then gets the options for a driver, parses them from
    command line and runs the driver"""
    parser = populate_parser(exec_drivers)
    args = parser.parse_args(sysargs[1:])
    args.func(**vars(args))
