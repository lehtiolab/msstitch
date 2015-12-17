def start_msstitch(drivers, passed_command, args):
    commandmap = {driver.get_commandname() for driver in drivers}
    command = commandmap[passed_command](**vars(args))
    command.run()
