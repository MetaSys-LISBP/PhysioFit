"""
Module containing the Command-Line Interface control logic.

The Command-Line Interface has two goals:

    i) Have a local CLI that can be used to launch jobs in a concise and quick
    manner, without use of the Graphical User Interface. This can help when
    testing many different parameter sets, and also permits easy pipelining
    of the tool, as described below.
    ii) Have a version compatible with an integration onto the Galaxy W4M
    platform to be used in automated flux calculation workflows

"""
import argparse
import logging
from pathlib import Path
import sys

import pandas as pd

from physiofit.base.io import IoHandler, StandardDevs, ConfigParser

logger = logging.getLogger("physiofit")
logger.setLevel(logging.DEBUG)


def args_parse():
    """
    Parse arguments from user input.
    :return: Argument Parser object
    :rtype: class: argparse.ArgumentParser
    """

    parser = argparse.ArgumentParser(
        "Physiofit: Extracellular flux estimation software"
    )

    # Parse data arguments (tsv + json)
    parser.add_argument(
        "-d", "--data", type=str,
        help="Path to data file in tabulated format (txt or tsv)"
    )
    parser.add_argument(
        "-c", "--config", type=str,
        help="Path to config file in yaml format"
    )
    parser.add_argument(
        "-m", "--model", type=str,
        help="Which model should be chosen. Useful only if generating "
             "related config file"
    )

    # Parse developer arguments
    parser.add_argument(
        "-g", "--galaxy", action="store_true",
        help="Is the CLI being used on the galaxy platform"
    )
    parser.add_argument(
        "--list", action="store_true",
        help="Return the list of models in model folder"
    )
    parser.add_argument(
        "-v", "--debug_mode", action="store_true",
        help="Activate the debug logs"
    )

    # Parse selective output path arguments (for galaxy implementation mostly)
    # parser.add_argument(
    #     "-op", "--output_pdf", type=str,
    #     help="Path to output the pdf file containing plots"
    # )
    # parser.add_argument(
    #     "-of", "--output_fluxes", type=str,
    #     help="Path to output the flux results"
    # )
    # parser.add_argument(
    #     "-os", "--output_stats", type=str,
    #     help="Path to output the khi² test"
    # )
    parser.add_argument(
        "-oc", "--output_config", type=str,
        help="Path to output the yaml config file"
    )
    parser.add_argument(
        "-or", "--output_recap", type=str,
        help="Path to output the summary"
    )
    parser.add_argument(
        "-oz", "--output_zip", type=str,
        help="Path to export zip file"
    )

    return parser


def run(data, args, experiments):
    """
    Run the optimization process
    :param data: Dataframe containing the main dataset
    :type data: pandas.DataFrame
    :param args: Arguments from the user input
    :type args: argparse.Namespace
    :param experiments: List of experiments to process
    :type experiments: list
    """

    for exp in experiments:
        io = IoHandler()
        logger.info(f"Processing experiment: {exp}")
        exp_data = data.loc[exp, :].sort_values("time")
        exp_data = exp_data.reset_index().drop("experiments", axis=1).copy()
        logger.info(f"Input Data: \n{exp_data}")
        if args.galaxy:
            io.wkdir = Path('.')
        else:
            io.wkdir = Path(args.data).parent
        logger.info(f"Home path: {io.wkdir}")
        logger.info(f"Reading configuration file at {args.config}")
        io.configparser = io.read_yaml(args.config)
        logger.info(f"Config File:\n{io.configparser}")
        model = io.select_model(io.configparser.model["model_name"], exp_data)
        model.get_params()
        logger.info(f"Selected Model:\n{model}")
        model = io.configparser.update_model(model)
        logger.info(f"Updated Model: \n{model}")
        logger.info(f"Model Data: \n{model.data}")

        fitter = io.initialize_fitter(
            model.data,
            model=model,
            mc=io.configparser.mc,
            iterations=io.configparser.iterations,
            sd=io.configparser.sds,
            debug_mode=args.debug_mode
        )
        fitter.optimize()
        if fitter.mc:
            fitter.monte_carlo_analysis()
        fitter.khi2_test()
        try:
            fitter.aic_test()
        except ValueError:
            logger.warning(
                "Not enough measurements to calculate AIC"
            )
            fitter.aic, fitter.aic_c = "NA"
        df = pd.DataFrame.from_dict(
            fitter.parameter_stats,
            orient="columns"
        )
        df.index = [
            f"{exp} {param}" for param in fitter.model.parameters.keys()
        ]
        if not io.multiple_experiments:
            io.multiple_experiments = []
        io.multiple_experiments.append(df)
        res_path = io.wkdir / (io.wkdir.name + "_res") / exp
        logger.info(res_path)
        if not res_path.is_dir():
            res_path.mkdir(parents=True)
        logger.info(f"Results:\n{df}")
        io.output_report(fitter, str(res_path))
        io.output_plots(fitter, str(res_path))
        io.output_pdf(fitter, str(res_path))
        io.figures = []
    if args.output_zip:
        output_dir = res_path.parents[0]
        generate_zips(str(output_dir), args.output_zip)
    logger.debug(f"Dataframes to concatenate:\n{io.multiple_experiments}")

    if args.galaxy:
        io.output_recap(export_path=args.output_recap, galaxy=args.galaxy)
    else:
        io.output_recap(export_path=str(res_path.parent), galaxy=False)


def generate_zips(path_to_data_folder, output_path):
    """Generate output zip file containing results

    :param path_to_data_folder: Path to the data folder
    :type path_to_data_folder: str
    :param output_path: Path to the output zip file
    :type output_path: str
    """
    from zipfile import ZipFile
    directory = Path(path_to_data_folder)
    output = Path(output_path)
    with ZipFile(str(output), mode="w") as archive:
        for file_path in directory.rglob("*"):
            archive.write(
                file_path,
                arcname=file_path.relative_to(directory)
            )
    return


def generate_config(args, data, _logger):
    """"
    Generate a configuration file for the selected model
    :param args: Arguments from the user input
    :type args: argparse.Namespace
    :param data: Dataframe containing the main dataset
    :type data: pandas.DataFrame
    :param _logger: Logger object
    :type _logger: logging.Logger
    """
    _logger.info("Launching in configuration file generator mode")
    # Run checks
    if args.model is None:
        # Must be path to model.py file or name of a model
        raise ValueError(
            "Please select a model to generate the associated configuration "
            "file"
        )
    if args.output_config is None:
        raise ValueError(
            "Please select an output directory for the configuration file"
        )
    try:
        io = IoHandler()
        if Path(args.model).is_file():
            model_class = io.read_model(args.model)
            model = model_class(data)
            model.get_params()

        else:
            model = io.select_model(args.model, data)
            model.get_params()
    except Exception:
        _logger.error("There was an error while initialising the model")
        raise

    # Build config parser and export configuration file
    default_sds = StandardDevs({col: 0.2 for col in model.name_vector})
    configparser = ConfigParser(
        path_to_data=args.data,
        selected_model=model,
        sds=default_sds,
        mc=True,
        iterations=100
    )
    configparser.export_config(args.output_config)
    _logger.info(
        f"Finished exporting the configuration file to {args.output_config}")
    sys.exit()


def _build_logger(mode):
    """
    Build the logger object
    :param mode: Debug mode
    :type mode: bool
    """
    _logger = logging.getLogger("physiofit")
    cli_handle = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - '
                                  '%(message)s')
    cli_handle.setFormatter(formatter)
    if mode:
        cli_handle.setLevel(logging.DEBUG)
    else:
        cli_handle.setLevel(logging.INFO)
    _logger.addHandler(cli_handle)
    return _logger


def process(args):
    """
    Process the arguments
    :param args: Arguments from the user input
    :type args: argparse.Namespace
    :return: None
    """
    _logger = _build_logger(args.debug_mode)
    _logger.info("Starting Physiofit CLI")
    if args.list:
        IoHandler.get_model_list()
        sys.exit()

    if args.data is None:
        raise ValueError("A path towards the data must be given")

    # Ensure that the path to data is correct
    path_to_data = Path(args.data)
    if not path_to_data.is_file():
        raise ValueError(
            f"The data path is not correct. Please check and try again. "
            f"Input data path: \n{args.data}"
        )
    # Ensure that the input file is a tsv if we are local
    if not args.galaxy:
        if path_to_data.suffix not in [".tsv", ".txt"]:
            raise TypeError(
                f"The input data must be in tsv/txt format. Detected format:"
                f" {path_to_data.suffix}"
            )
    else:
        if not path_to_data.suffix == ".dat":
            raise TypeError(
                f"The input data must be in dat format (galaxy data "
                f"format). Detected format: {path_to_data.suffix}"
            )

    # Read & check data
    data = IoHandler.read_data(str(path_to_data))

    # If no configuration file is present we assume the user wants to
    # generate a default one
    if args.config is None:
        generate_config(args, data, _logger)

    # If configuration file is present we launch the optimization
    experiments = list(data["experiments"].unique())
    data = data.set_index("experiments")
    run(data, args, experiments)
    _logger.info("Done!")
    sys.exit()


def main():
    """Main routine for the CLI interface"""
    parser = args_parse()
    args = parser.parse_args()
    process(args)
