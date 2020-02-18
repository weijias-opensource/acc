# -*- coding: utf-8 -*-
import click

from acc.io import import_data
from acc.processing import processing_event_Z
from acc.processing import processing_noise_Z
from acc.migration import migration_1station
from acc.plotting import plot_profile
from acc.util import pkl2sac1

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS)
def run():
    pass


@run.command()
@click.option("-p", "--paras", default="conf.json", help="The parameter JSON file")
def importdata(**kwargs):
    """Import data and convert to internal format."""
    import_data(kwargs["paras"])


@run.command()
@click.option("-p", "--paras", default="conf.json", help="The parameter JSON file")
def calevent(**kwargs):
    """Calculate teleseismic P-wave autocorrelograms."""
    processing_event_Z(kwargs["paras"])


@run.command()
@click.option("-p", "--paras", default="conf.json", help="The parameter JSON file")
def calnoise(**kwargs):
    """Calculate ambient noise autocorrelograms."""
    processing_noise_Z(kwargs["paras"])


@run.command()
@click.option("-p", "--paras", default="conf.json", help="The parameter JSON file")
def migration(**kwargs):
    """Migrate teleseismic P-wave autocorrelograms."""
    migration_1station(kwargs["paras"])


@run.command()
@click.option("-p", "--paras", default="conf.json", help="The parameter JSON file")
def plotprofile(**kwargs):
    """Plot migrated P-wave profile."""
    plot_profile(kwargs["paras"])


@run.command()
@click.option("-d", "--directory", help="The direcotry containing files to be converted")
@click.option("-f", "--format", default="SAC", help="The target format, support SAC, MSEED. default: SAC")
@click.option("-s", "--suffix", default="pkl", help="The suffix of files to be converted. default: pkl")
def pkl2sac(**kwargs):
    """Convert pkl format to SAC or MSEED format."""
    print(kwargs)
    pkl2sac1(directory=kwargs["directory"], suffix=kwargs["suffix"], fmt=kwargs["format"])


if __name__ == '__main__':
    run()
