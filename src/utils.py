import logging
def get_logger(__name__):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)
    try:
        fh = logging.FileHandler(f"{__file__}.log")
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    except PermissionError:
        pass

    # fh.setLevel(logger.DEBUG)
    # create console handler with a higher log level
    return logger


logger = get_logger(__name__)