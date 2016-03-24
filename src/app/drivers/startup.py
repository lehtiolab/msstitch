import os
from argparse import ArgumentParser, RawTextHelpFormatter

VERSION_NUMBER = '1.0'


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


def set_command_on_parser(drivers, commandmap, parser):
    commandhelptxt = ['\033[1m{}\033[0m\n{}'.format(name,
                                                    driver.get_commandhelp())
                      for name, driver in commandmap.items()]
    parser.add_argument(dest='command', type=str,
                        help='How to manipulate the input:'
                        '\n\n{}'.format('\n\n'.join(commandhelptxt)))
    return commandmap


def get_cmd_driver(drivers, commandline):
    cmdparser = ArgumentParser(formatter_class=RawTextHelpFormatter)
    commandmap = {driver.get_commandname(): driver for driver in drivers}
    set_command_on_parser(drivers, commandmap, cmdparser)
    try:
        passed_command = [commandline[1]]
    except IndexError:
        passed_command = []
    cmdargs = cmdparser.parse_args(passed_command)
    try:
        driver = commandmap[cmdargs.command]
    except KeyError:
        # Can we use argparse to deliver a friendly message?
        raise RuntimeError('Command not recognized: {}'.format(
            cmdargs.command))
    else:
        return driver, passed_command


def start_msstitch(exec_drivers, sysargs):
    """Passed all drivers of executable, checks which command is passed to
    the executable and then gets the options for a driver, parses them from
    command line and runs the driver"""
    parser = populate_parser(exec_drivers)
    args = parser.parse_args(sysargs[1:])
    args.func(**vars(args))
