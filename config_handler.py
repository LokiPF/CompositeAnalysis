import configparser

from configuration import Configuration


def read_config(config_dataclass: Configuration) -> bool:
    """
    Reads the config file and returns True if the config file is valid.
    :return: bool
    """
    try:
        print("Reading config file...", end='')
        config = configparser.ConfigParser()
        config.read('config.ini')
    except configparser.Error:
        print("Error!")
        return False
    print("OK!")
    set_values(config=config, cd=config_dataclass)
    return True


def set_values(config: configparser.ConfigParser(), cd: Configuration) -> bool:
    """
    Takes a config file and sets the values in the dataclasses.
    :return: bool
    """
    print("Setting Dataclasses...", end='')
    try:
        cd.knockdown_factor = float(config['configuration']['knockdown_factor'])
        cd.sigma_ul = float(config['configuration']['sigma_ul'])
        cd.ul_factor = float(config['configuration']['ul_factor'])
    except KeyError:
        print("Error!")
        return False
    print("OK!")
    return True
