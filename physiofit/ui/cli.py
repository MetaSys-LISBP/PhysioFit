"""
Module containing the Command-Line Interface control logic.

The Command-Line Interface has two goals:

    i) Have a local CLI that can be used to launch jobs in a concise and quick manner, without use of the Graphical
    User Interface. This can help when testing many different parameter sets, and also permits easy pipelining of the
    tool, as described below.
    ii) Have a version compatible with an integration onto the Galaxy W4M platform to be used in automated flux
    calculation workflows

"""
import argparse
import logging
from pathlib import Path
import sys
# import shutil

import pandas as pd

from physiofit.base.io import IoHandler, StandardDevs, ConfigParser

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
        help="Which model should be chosen. Useful only if generating related config file"
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
    parser.add_argument(
        "-op", "--output_pdf", type=str,
        help="Path to output the pdf file containing plots"
    )
    parser.add_argument(
        "-of", "--output_fluxes", type=str,
        help="Path to output the flux results"
    )
    parser.add_argument(
        "-os", "--output_stats", type=str,
        help="Path to output the khi² test"
    )
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

def run(data, args, logger, experiments):

    for exp in experiments:
        io = IoHandler()
        exp_data = data.loc[exp, :].sort_values("time").copy()
        logger.info(f"Input Data: \n{exp_data}")
        io.home_path = Path(args.data).parent
        logger.info(f"Home path: {io.home_path}")
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
        df = pd.DataFrame.from_dict(
            fitter.parameter_stats,
            orient="columns"
        )
        df.index = [
            f"{exp} {param}" for param in fitter.model.parameters_to_estimate.keys()
        ]
        if not io.multiple_experiments:
            io.multiple_experiments = []
        io.multiple_experiments.append(df)
        if exp is not None:
            res_path = io.home_path / (io.home_path.name + "_res") / exp
        else:
            res_path = io.home_path / (io.home_path.name + "_res")
        logger.info(res_path)
        if not res_path.is_dir():
            res_path.mkdir(parents=True)
        logger.info(f"Results:\n{df}")
        io.output_report(fitter, str(res_path))
        io.output_plots(fitter, str(res_path))
        io.output_pdf(fitter, str(res_path))
        io.figures = []
    if args.output_zip:
        # logger.info(io.home_path)
        # logger.info(io.home_path / io.home_path.name / "_res")
        # logger.info(args.output_zip)
        dir = res_path.parents[0]
        # for line in dir.rglob("*"):
        #     print(line)
        generate_zips(str(dir), args.output_zip, logger)

    # shutil.make_archive(
    #     args.output_zip,
    #     format='zip',
    #     root_dir=str(io.home_path / (io.home_path.name + "_res"))
    #     )
    io.output_recap(export_path=args.output_recap, galaxy=True)

def generate_zips(path_to_data_folder, output_path, logger):
    """Generate output zip file containing results

    Args:
        path_to_data_folder (str): path to folder containing directories & files to zip
        output_path (str): path to export archive to
        logger (Logging.logger): main logger

    Returns:
        None:
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
    return 0

def generate_config(args, data, logger):

    logger.info("Launching in configuration file generator mode")
    # Run checks
    if args.model is None:
        # Must be path to model.py file or name of a model
        raise ValueError(
            f"Please select a model to generate the associated configuration file"
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
        logger.error("There was an error while initialising the model")
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
    logger.info(f"Finished exporting the configuration file to {args.output_config}")
    sys.exit()

def process(args):

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    out_stream = logging.StreamHandler()
    formatter = logging.Formatter()
    out_stream.setFormatter(formatter)
    if args.debug_mode:
        out_stream.setLevel(logging.DEBUG)
    else:
        out_stream.setLevel(logging.INFO)
    logger.addHandler(out_stream)

    if args.list:
        IoHandler.get_model_list()
        sys.exit()

    if args.data is None:
        raise ValueError("A path towards the data must be given")

    # Ensure that the path to data is correct
    path_to_data = Path(args.data)
    if not path_to_data.is_file():
        raise ValueError(
            f"The data path is not correct. Please check and try again. Input data path: \n{args.data}"
        )
    # Ensure that the input file is a tsv if we are local
    if not args.galaxy:
        if not path_to_data.suffix in [".tsv", ".txt"]:
            raise TypeError(
                f"The input data must be in tsv/txt format. Detected format: {path_to_data.suffix}"
            )
    else:
        if not path_to_data.suffix == ".dat":
            raise TypeError(
                f"The input data must be in dat format (galaxy data format). Detected format: {path_to_data.suffix}"
            )

    # Read & check data
    data = IoHandler.read_data(str(path_to_data))

    # If no configuration file is present we assume the user wants to generate a default one
    if args.config is None:
        generate_config(args, data, logger)

    # If configuration file is present we launch the optimization
    experiments = list(data["experiments"].unique())
    data = data.set_index("experiments")
    for exp in experiments:
        logger.info(f"Running optimization for {exp}")
        run_data = data.loc[exp, :].sort_values("time").copy()
        run(run_data, args, logger, exp)
    run(data, args, logger, experiments)
    logger.info("Done!")
    sys.exit()

def main():
    parser = args_parse()
    args = parser.parse_args()
    process(args)
