def start_msstitch(drivers, parser, sysargs):
    passed_command = sysargs[1]
    args = parser.parse_args(sysargs[2:])
    commandmap = {driver.get_commandname() for driver in drivers}
    command = commandmap[passed_command](**vars(args))
    command.run()
